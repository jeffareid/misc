/*----------------------------------------------------------------------*/
/*      eccr16z.c       ecc rs GF(2^4) + syndromes + handle zero        */
/*                                                                      */
/*      Jeff Reid       2025AUG15 15:00                                 */
/*----------------------------------------------------------------------*/
#define _CRT_SECURE_NO_WARNINGS 1       /* disable sscanf warnings */

#include <stdio.h>
#include <string.h>
#include <memory.h>

/* only set one of these to 1 */
/* Berelekamp Welch decode */
#define BW 0
/* Gao decode */
#define GAO 0
/* syndrome decode */
#define SYN 1

/* powers of alpha from 0 to n-2, not used */
#define FOURIER 0

/* allow user to input polynomial, not used */
#define ENCPOLY 0

typedef unsigned char BYTE;

#define MAXDAT 16
#define MAXPAR 12

/* display matrix inv stuff */
#define DISPLAYI 0
/* display vandermonde stuff */
#define DISPLAYV 0
/* display matrix stuff */
#define DISPLAYM 0
/* display euclid stuff */
#define DISPLAYE 0
/* display syndrome stuff */
#define DISPLAYS 1
#define DISPLAYX 1
/* display decode stuff */
#define DISPLAYD 0

typedef struct{                         /* vector structure */
    int  size;
    BYTE data[18];                       /*  +2 for lagrange */
}VECTOR;
/* POLY == VECTOR */
#define POLY VECTOR

typedef struct{                         /* euclid structure */
    int  size;                          /* # of data values */
    int  indx;                          /* index to right side */
    BYTE data[20];                      /* left and right side data */
}EUCLID;

typedef struct{                         /* matrix structure */
    int  nrows;
    int  ncols;
    BYTE data[16*16*2];
}MATRIX;

/*----------------------------------------------------------------------*/
/*      data                                                            */
/*----------------------------------------------------------------------*/
static int aiGA[6] =                    /* GF's and Alpha's */
   {0x13, 0x02, 0x19, 0x02, 0x1f, 0x03};

static int      iGF;                    /* GF(2^8) poly */
static BYTE     bAlpha;
static BYTE     bPad0[3];

static BYTE     abExp[32];
static BYTE     abLog[16];
static BYTE     abInv[16];
static BYTE     abNeg[16];
static int      f[16];                  /* fix indexes */

static BYTE     abId[32];               /* array for MatrixInv */

static VECTOR   vData;                  /* data */
static VECTOR   vPoly;                  /* poly */

static VECTOR   vLagrK;                 /* lagrange stuff */
static VECTOR   vLagrN;
static VECTOR   vLagrs;
static VECTOR   vLagr;

#if FOURIER
static VECTOR   vTran;                  /* transform stuff */
static VECTOR   vITrn;
#endif

static VECTOR   vA;                     /* evaluation points */
static VECTOR   vA2I;                   /* A to index of A */

#if SYN
/*                                      ** for syndrome based correction */
static VECTOR   vB;                     /* subset of A */
static VECTOR   vL;                     /* evaluations of L[] */
static VECTOR   vU;                     /* 1/vL[]() */
static VECTOR   vS;                     /* syndromes */
static VECTOR   pE;                     /* error locator polynomial */
static VECTOR   pD;                     /* derivative of pE */
static VECTOR   vI;                     /* inverse locators */
static VECTOR   vO;                     /* offsets */
static VECTOR   vX;                     /* vA values corresponding to offsets */
static VECTOR   vY;                     /* vU values corresponding to offsets */
static VECTOR   vV;                     /* error values */
static VECTOR   pO;                     /* pOmega */
static VECTOR   vF;                     /* error values */
static VECTOR   pP[16];                 /* polynomials */
static VECTOR   pL[16];                 /* L[] polynomials */
static VECTOR   pS[16];                 /* S[] syndrome polynomials */
#endif

static VECTOR   vE;                     /* error poly */
static VECTOR   vP;                     /* mapping poly */
static VECTOR   vQ;                     /* vP*vE */
static MATRIX   mM;                     /* matrix */

static MATRIX   mV;                     /* Vandermonde matrix */
static MATRIX   mK;                     /* Vandermode matrix k by k */
static MATRIX   mI;                     /* inverse of Mk */
static MATRIX   mS;                     /* Vandermonde matrix systematic */
static MATRIX   mD;                     /* data */
static MATRIX   mE;                     /* encoded data */

static EUCLID   E0;                     /* for DecodeG() and DecodeS() */
static EUCLID   E1;

static int      E;                      /* max number of errors */
static int      F;                      /* N - K */
static int      K;                      /* number data - parities */
static int      N;                      /* number data + num par */

static BYTE     abUser[80];             /* user input buffer */
/*----------------------------------------------------------------------*/
/*      code                                                            */
/*----------------------------------------------------------------------*/
static void     LagrangeI(VECTOR *, int);
static void     Lagrange(VECTOR *, int);
#if FOURIER
static void     LagrangeIP(void);
static void     Fourier(void);
static void     FourierI(void);
#endif
static void     Encode(void);
static void     EncodeP(void);
static void     EncodeS(void);
static void     EncodeV(void);
static void     DecodeBW(void);
static void     DecodeG(void);
static void     DecodeS(void);
static void     Fix(void);
static void     PolyDiv(POLY *, POLY *, POLY *);
static void     MatrixMpy(MATRIX *, MATRIX *, MATRIX *);
static int      MatrixInv(MATRIX *, MATRIX *);
static int      Poly2Root(VECTOR *, VECTOR *);
static void     Root2Poly(VECTOR *, VECTOR *);
static BYTE     Poly2Val(VECTOR *, BYTE);
static BYTE     GFAdd(BYTE, BYTE);
static BYTE     GFSub(BYTE, BYTE);
static BYTE     GFMpy(BYTE, BYTE);
static BYTE     GFDiv(BYTE, BYTE);
static BYTE     GFPow(BYTE, int);
static BYTE     GFMpy0(BYTE, BYTE);
static void     ShowVector(VECTOR *);
static void     ShowEuclid(EUCLID *);
static void     ShowMatrix(MATRIX *);
static int      Conrs(void);
static void     DoUser(void);

/*----------------------------------------------------------------------*/
/*      Vandermonde  generate Vandermonde matrices                      */
/*----------------------------------------------------------------------*/
static void Vandermonde()
{
int r, c;
    mV.nrows = K;
    mV.ncols = N;
    for(r = 0; r < K; r++)
        for(c = 0; c < N; c++)
            mV.data[r*N+c] = GFPow(vA.data[c], r);
#if DISPLAYV
    printf("Vandermonde\n");
    ShowMatrix(&mV);
#endif
    mK.nrows = K;
    mK.ncols = K;
    for(r = 0; r < K; r++)
        for(c = 0; c < K; c++)
            mK.data[r*K+c] = mV.data[r*N+c];
#if DISPLAYV
    printf("Vandermonde k x k\n");
    ShowMatrix(&mK);
#endif
    if(MatrixInv(&mI, &mK)){
        printf("can't invert\n");
        return;
    }
#if DISPLAYV
    printf("inverse\n");
    ShowMatrix(&mI);
#endif
    MatrixMpy(&mS, &mI, &mK);
#if DISPLAYV
    printf("check\n");
    ShowMatrix(&mS);
#endif
    MatrixMpy(&mS, &mI, &mV);
#if DISPLAYV
    printf("systematic\n");
    ShowMatrix(&mS);
#endif
}

/*----------------------------------------------------------------------*/
/*      LagrangeI   initialize Lagrange dvnd                            */
/*----------------------------------------------------------------------*/
static void LagrangeI(VECTOR *pL, int S)
{
int i, j, ai;

    S += 1;                         /* pL = II(for i = 0 to S)(x-a(i)) */
    pL->size = S;
    for(i = 1; i < S; i++)          /* zero pL */
        pL->data[i] = 0;
    pL->data[0] = 1;                /* ms coef pL == 1 */
    ai = 0;                         /* index for a() */
    for(j = 1; j < S; j++){         /* generate pL numerator */
        for(i = j; i > 0; i--)
            pL->data[i] = GFSub(pL->data[i],
                GFMpy(pL->data[i-1], vA.data[ai]));
        ai++;
    }
}

/*----------------------------------------------------------------------*/
/*      Lagrange   interpolate S term polynomial                        */
/*----------------------------------------------------------------------*/
static void Lagrange(VECTOR *pL, int S)
{
int i, j;
BYTE t;
    vLagrs.size = S;
    vLagr.size = S;
    for(i = 0; i < S; i++)              /* reset vLagr */
        vLagr.data[i] = 0;
    for(j = 0; j < S; j++){
        vLagrs.data[0] = t = 1;
        for(i = 1; i < S; i++){         /* vLagrs = pL/(1 x - a(j)) */
            t = GFAdd(pL->data[i], GFMpy(t,vA.data[j]));
            vLagrs.data[i] = t;
        }
        t = 1;                          /* generate denominator */
        for(i = 0; i < S; i++){
            if(i == j)
                continue;
            t = GFMpy(t, GFSub(vA.data[j], vA.data[i]));   /* check this */
        }
        for(i = 0; i < S; i++)          /* (vLagrs / t) * vData[j] */
            vLagrs.data[i] = GFMpy(GFDiv(vLagrs.data[i], t), vData.data[j]);
        for(i = 0; i < S; i++)          /* vLagr += vLagrs */
            vLagr.data[i] = GFAdd(vLagr.data[i], vLagrs.data[i]);
    }
}

/*----------------------------------------------------------------------*/
/*      Lagrange interpolate values to coefficients                     */
/*          doesn't require sequential powers of alpha                  */
/*----------------------------------------------------------------------*/
static void LagrangeIP(void)
{
int i;

    Lagrange(&vLagrN, vData.size);
    for(i = 0; i < vData.size; i++)
        vData.data[i] = vLagr.data[i];
}

#if FOURIER
/*----------------------------------------------------------------------*/
/*      Fourier     Fourier transform coefficients to values            */
/*                  same calculation as syndrome                        */
/*----------------------------------------------------------------------*/
static void Fourier(void)
{
int i,j;

    vTran.size = N;                     /* set size */
    for(j = 0; j < vTran.size; j++){
        vTran.data[j] = vData.data[0];  /* generate a transform */
        for(i = 1; i < vData.size;  i++){
            vTran.data[j] = GFAdd(vData.data[i],
                GFMpy(a(j), vTran.data[j]));}}
    for(i = 0; i < vData.size; i++)
        vData.data[i] = vTran.data[i];
}

/*----------------------------------------------------------------------*/
/*      FourerI     inverse Fourier transform                           */
/*----------------------------------------------------------------------*/
static void FourierI(void)
{
int i, j;
BYTE t;

    vITrn.size = N;                     /* encode data to vITrn */
    for(j = 0; j < N; j++){
        t = vData.data[N-1];
        for(i = N-2; i >=0; i--)
            t = GFAdd(GFMpy(t, a(N-j)), vData.data[i]);
#if 1
        vITrn.data[j] = t;              /* N mod 2 == 1*/
#else
        vITrn.data[j] = GFDiv(t, N);    /* N mod p for GF(p^1) */
#endif
    }
    for(j = 0; j < N; j++)
        vData.data[N-j-1] = vITrn.data[j];
}
#endif

/*----------------------------------------------------------------------*/
/*      Encode                                                          */
/*----------------------------------------------------------------------*/
static void Encode(void)
{
int i, j;
BYTE t;

    Lagrange(&vLagrK, K);               /* interpolate to vLagr */
    vPoly.size = vLagr.size;            /* copy to vPoly */
    for(i = 0; i < vLagr.size; i++)
        vPoly.data[i] = vLagr.data[i];

    for(j = 0; j < N; j++){             /* check result | encode */
        t = vPoly.data[0];
        for(i = 1; i < K; i++)
            t = GFAdd(GFMpy(t, vA.data[j]), vPoly.data[i]);
        if(j < K && t != vData.data[j]){
            printf("encode error\n");
            break;
        }
        vData.data[j] = t;
    }
}

/*----------------------------------------------------------------------*/
/*      EncodeS                                                         */
/*----------------------------------------------------------------------*/
static void EncodeS(void)
{
int i;
    mD.nrows = 1;
    mD.ncols = K;
    for(i = 0; i < K; i++)
        mD.data[i] = vData.data[i];
    MatrixMpy(&mE, &mD, &mS);
    for(i = 0; i < N; i++)
        vData.data[i] = mE.data[i];
}

/*----------------------------------------------------------------------*/
/*      EncodeV                                                         */
/*----------------------------------------------------------------------*/
static void EncodeV(void)
{
int i;
    mD.nrows = 1;
    mD.ncols = K;
    for(i = 0; i < K; i++)
        mD.data[i] = vData.data[i];
    MatrixMpy(&mE, &mD, &mV);
    for(i = 0; i < N; i++)
        vData.data[i] = mE.data[i];
}

#if ENCPOLY
/*----------------------------------------------------------------------*/
/*      EncodeP                                                         */
/*----------------------------------------------------------------------*/
static void EncodeP(void)
{
int i, j;
BYTE t;

    for(j = 0; j < N; j++){             /* encode from poly */
        t = vPoly.data[0];
        for(i = 1; i < K; i++)
            t = GFAdd(GFMpy(t, vA.data[j]), vPoly.data[i]);  /* check this */
        vData.data[j] = t;
    }
}
#endif

/*----------------------------------------------------------------------*/
/*      UnEncodeV                                                       */
/*----------------------------------------------------------------------*/
static void UnEncodeV(void)
{
int i;
    Lagrange(&vLagrN, vData.size);
    for(i = 0; i < vData.size; i++)
        vData.data[i] = vLagr.data[N-i-1];
}

/*----------------------------------------------------------------------*/
/*      DecodeBW    Berlekamp Welch  P() = Q()/E(()                     */
/*----------------------------------------------------------------------*/
static void DecodeBW(void)
{
int i, j, k, e;
BYTE *prowj, *prowk;
BYTE t;
    e = E;                              /* assume max errors */
    mM.nrows = N;
    mM.ncols = N+1;
loop:
    /* generate matrix */
    for(j = 0; j < mM.nrows; j++){
        prowj = mM.data+j*mM.ncols;
        for(i = 0; i < e; i++)
            prowj[i] = GFMpy(vData.data[j],GFPow(vA.data[j],i));  /* check */
        for(; i < mM.ncols-1; i++)
            prowj[i] = GFSub(0, GFPow(vA.data[j], (i-e)));
        prowj[i] = GFSub(0, GFMpy(vData.data[j], GFPow(vA.data[j],e)));
    }
#if DISPLAYM
    printf("\n"); ShowMatrix(&mM);
#endif
    /* solve matrix */
    for(k = 0; k < mM.nrows; k++){
        prowk = mM.data+k*mM.ncols;
        for(j = k; j < mM.nrows; j++)   /* find non-zero in column */
            if(mM.data[j*mM.ncols+k] != 0)
                break;
        if(j == mM.nrows){              /* if redundant try 1 less error */
            if(--e < 0)                 /*  if 0 errors return */
                return;
#if 0                                   /* next two lines are optional */
            mM.nrows -= 2;              /*  Q() will just end up with  */
            mM.ncols -= 2;              /*  leading zeroes with #if 0   */
#endif
            goto loop;
        }
        if(j != k){                     /* swap rows if needed */
            prowj = mM.data+j*mM.ncols;
            for(i = 0; i < mM.ncols; i++){
                t        = prowk[i];
                prowk[i] = prowj[i];
                prowj[i] = t;
            }
        }
        t = prowk[k];                   /* divide by M[k][k] */
        for(i = 0; i < mM.ncols; i++)
            prowk[i] = GFDiv(prowk[i], t);
        for(j = 0; j < mM.nrows; j++){  /* zero out columns */
            if(j == k)
                continue;
            prowj = mM.data+j*mM.ncols;
            t = prowj[k];
            for(i = 0; i < mM.ncols; i++)
                prowj[i] = GFSub(prowj[i], GFMpy(prowk[i], t));
        }
    }
#if DISPLAYM
    printf("\n"); ShowMatrix(&mM);
#endif
    vE.size = e+1;
    vQ.size = mM.nrows-e;
    /* set proj to 1 row past last int of matrix */
    prowj = mM.data+((mM.nrows+1)*mM.ncols-1);
    for(i = 0; i < vQ.size; i++){       /* create vQ */
        prowj -= mM.ncols;
        vQ.data[i] = *prowj;
    }
    vE.data[0] = 1;                     /* create vE */
    for(i = 1; i < vE.size; i++){
        prowj -= mM.ncols;
        vE.data[i] = *prowj;
    }
#if DISPLAYM
    printf("\n");
#endif
#if DISPLAYD
    printf("vQ:\n"); ShowVector(&vQ);
#endif
    while(vQ.size != 0 && vQ.data[0] == 0){   /* remove leading zeroes */
        vQ.size--;
        memcpy(vQ.data, vQ.data+1, vQ.size*sizeof(BYTE));
    }
#if DISPLAYD
    printf("vQ:\n"); ShowVector(&vQ);
    printf("vE:\n"); ShowVector(&vE);
#endif
    PolyDiv(&vP, &vQ, &vE);
#if DISPLAYD
    printf("vQ:\n"); ShowVector(&vQ);
    printf("vP:\n"); ShowVector(&vP);
#endif
}

/*----------------------------------------------------------------------*/
/*      DecodeG     Gao (Euclid) decoder                                */
/*----------------------------------------------------------------------*/
static void DecodeG(void)
{
EUCLID  *pE0;                           /* R[i-2] / A[i-1] */
EUCLID  *pE1;                           /* R[i-1] / A[i-2] */
EUCLID  *pET;                           /* temp */
int     i, j;
BYTE    t;

/*      E0 = R[0] | A[1]                    */
/*      R[0] = product (x-a(i)) | A[1] = 1  */
    E0.size = N+2;
    E0.indx = N+1;
    for(i = 1; i < E0.indx; i++)        /* zero R[0] */
        E0.data[i] = 0;
    E0.data[0] = 1;                     /* R[0][0] = 1 */
    for(j = 1; j < E0.indx; j++){       /* generate rest of R[0] */
        for(i = j; i > 0; i--)
            E0.data[i] = GFSub(E0.data[i],
                GFMpy(E0.data[i-1], vA.data[j-1]));
    }
    E0.data[E0.indx] = 1;               /* A[1] = 1 */
    pE0 = &E0;

/*      E1 = R[1] | A[0]                    */
/*      R[1] = Lagrange polyn \ A[0] = 0    */

    Lagrange(&vLagrN, N);               /* interpolate to vLagr */
    E1.size = N+2;                      /* copy to E1 */
    E1.indx = N+1;
    E1.data[0] = 0;
    j = E1.indx - vLagr.size;
    for(i = 0; i < vLagr.size; i++)
        E1.data[j++] = vLagr.data[i];
    E1.data[E1.indx] = 0;               /* A[0] = 0 */
    pE1 = &E1;

/*      do Euclid algorithm */

    while(1){                           /* while deg of pE1[] >= (N+K)/2) */
#if DISPLAYE
        printf("E0: ");
        ShowEuclid(pE0);
        printf("E1: ");
        ShowEuclid(pE1);
#endif
        while((pE1->data[0] == 0) &&    /* shift dvsr left until msb!=0 */
              (pE1->indx != 0)){        /*  or fully shifted left */
            pE1->indx--;
            memcpy(&pE1->data[0], &pE1->data[1], (pE1->size-1)*sizeof(BYTE));
            pE1->data[pE1->size-1] = 0;
        }
        if(pE1->indx <= ((N+K)/2))      /* if (deg of pE1[] < (N+K)/2)break */
            break;
        while(1){                       /* while more divide sub-steps */
            if(pE0->data[0]){           /*   if dvnd msb!=0, do quot */
                t = GFDiv(pE0->data[0], pE1->data[0]);
                for(i = 0; i < pE1->indx; i++)
                    pE0->data[i] = GFSub(pE0->data[i], GFMpy(t, pE1->data[i]));
                for(i = pE0->indx; i < pE1->size; i++)
                    pE1->data[i] = GFSub(pE1->data[i], GFMpy(t, pE0->data[i]));
            }
            if(pE0->indx == pE1->indx)  /*   if divide done, break */
                break;
            pE0->indx--;                /*   shift dvnd */
            memcpy(&pE0->data[0], &pE0->data[1], (pE0->size-1)*sizeof(BYTE));
            pE0->data[pE0->size-1] = 0;
        }
        pET = pE1;                      /* swap pE0, pE1 */
        pE1 = pE0;
        pE0 = pET;
    }
#if DISPLAYE
    printf("X1: ");
    ShowEuclid(pE1);
#endif
    vQ.size = pE1->indx;                /* get vQ from pE1 */
    for(i = 0; i < vQ.size; i++)
        vQ.data[i] = pE1->data[i];
    vE.size = pE0->size - pE0->indx;    /* get vE from pE0 */
    j = pE0->size;
    for(i = 0; i < vE.size; i++)
        vE.data[i] = pE0->data[--j];
#if DISPLAYD
    printf("vQ:\n"); ShowVector(&vQ);
    printf("vE:\n"); ShowVector(&vE);
#endif
    t = vE.data[0];                     /* divide vE, vQ by vE[0] */
    if(t != 0){
        for(i = 0; i < vE.size; i++)    /*  to set vE[0] == 1 */
            vE.data[i] = GFDiv(vE.data[i], t);
        for(i = 0; i < vQ.size; i++)
            vQ.data[i] = GFDiv(vQ.data[i], t);
    }
#if DISPLAYD
    printf("vQ:\n"); ShowVector(&vQ);
    printf("vE:\n"); ShowVector(&vE);
#endif
    PolyDiv(&vP, &vQ, &vE);
#if DISPLAYD
    printf("vQ:\n"); ShowVector(&vQ);
    printf("vP:\n"); ShowVector(&vP);
#endif
}

#if SYN
/*----------------------------------------------------------------------*/
/*      DecodeS     syndrome decoder                                    */
/*----------------------------------------------------------------------*/
static void DecodeS(void)
{
EUCLID  *pE0;                           /* R[i-2] / A[i-1] */
EUCLID  *pE1;                           /* R[i-1] / A[i-2] */
EUCLID  *pET;                           /* temp */
int i, j;
BYTE    t;
    vS.size = F;                        /* generate syndromes */
    for(i = 0; i < F; i++)
        vS.data[i] = 0;
    for(j = 0; j < vData.size; j++){
        pS[j].size = F;
        for(i = 0; i < F; i++){
            pS[j].data[i] = GFMpy(pP[j].data[i], vData.data[j]);
            vS.data[F-1-i] = GFAdd(pS[j].data[i], vS.data[F-1-i]);
        }
    }
#if DISPLAYS
    printf("vS: ");
    ShowVector(&vS);
#endif
    t = 0;
    for(i = 0; i < vS.size; i++)
        if(vS.data[i] != 0)
            t = 1;
    if(t == 0){
        printf("no errors\n");
        return;
    }

/*      E0 = R[0] | A[1]                    */
/*      R[0] = z^F | A[1] = 1               */
    E0.size = F+2;
    E0.indx = F+1;
    for(i = 1; i < E0.indx; i++)        /* zero R[0] */
        E0.data[i] = 0;
    E0.data[0] = 1;                     /* R[0][0] = 1 */
    E0.data[E0.indx] = 1;               /* A[1] = 1 */
    pE0 = &E0;

/*      E1 = R[1] | A[0]                    */
/*      R[1] = S | A[0] = 0                 */

    E1.size = F+2;                      /* copy to E1 */
    E1.indx = F+1;
    E1.data[0] = 0;
    for(i = 0; i < F; i++)
        E1.data[i+1] = vS.data[i];
    E1.data[E1.indx] = 0;               /* A[0] = 0 */
    pE1 = &E1;

/*      do Euclid algorithm */

#if DISPLAYS
    printf("E0: ");
    ShowEuclid(pE0);
    printf("E1: ");
    ShowEuclid(pE1);
#endif
    while(1){                           /* while deg of pE1[] >= (N+K)/2) */
        while((pE1->data[0] == 0) &&    /* shift dvsr left until msb!=0 */
              (pE1->indx != 0)){        /*  or fully shifted left */
            pE1->indx--;
            memcpy(&pE1->data[0], &pE1->data[1], (pE1->size-1)*sizeof(BYTE));
            pE1->data[pE1->size-1] = 0;
        }
        if(pE1->indx <= E)              /* if (deg of pE1[] < E) break */
            break;
        while(1){                       /* while more divide sub-steps */
            if(pE0->data[0]){           /*   if dvnd msb!=0, do quot */
                t = GFDiv(pE0->data[0], pE1->data[0]);
                for(i = 0; i < pE1->indx; i++)
                    pE0->data[i] = GFSub(pE0->data[i], GFMpy(t, pE1->data[i]));
                for(i = pE0->indx; i < pE1->size; i++)
                    pE1->data[i] = GFSub(pE1->data[i], GFMpy(t, pE0->data[i]));
            }
            if(pE0->indx == pE1->indx)  /*   if divide done, break */
                break;
            pE0->indx--;                /*   shift dvnd */
            memcpy(&pE0->data[0], &pE0->data[1], (pE0->size-1)*sizeof(BYTE));
            pE0->data[pE0->size-1] = 0;
        }
        pET = pE1;                      /* swap pE0, pE1 */
        pE1 = pE0;
        pE0 = pET;
#if DISPLAYS
        printf("E0: ");
        ShowEuclid(pE0);
        printf("E1: ");
        ShowEuclid(pE1);
#endif
    }
    pE.size = pE0->size-pE0->indx;      /* set error poly size */
    for(i = 0; i < pE.size; i++)        /* pE = poly */
        pE.data[pE.size-1-i] = pE0->data[pE0->indx+i];
#if DISPLAYS
    printf("pE: ");
    ShowVector(&pE);
#endif
    t = pE.data[pE.size-1];             /* divide pE by lsc */
    if(t == 0){
        printf("least sig coef of pE == 0\n");
        return;}
    for(i = 0; i < pE.size; i++){
        pE.data[i] = GFDiv(pE.data[i], t);}
#if DISPLAYS
    printf("pE: ");
    ShowVector(&pE);
#endif

/*                                      ** derivative in GF(2^n) is just odd terms */
    pD.size = pE.size-1;                /* pD = derivative of pE */
    for(i = 0; i < pD.size; i++)
        pD.data[i] = GFMpy((BYTE)((pD.size-i)&1), pE.data[i]); /* GF(2^n), ...&1 */
#if DISPLAYS
    printf("pD: ");
    ShowVector(&pD);
#endif

    j = pE.size - 1;                    /* right shift E1 if Omega has leading zeroes */
    while(pE1->indx < j){
        pE1->indx++;
        for(i = pE1->size-1; i;){
            i--;
            pE1->data[i+1] = pE1->data[i];}
        pE1->data[0] = 0;}

    pO.size = pE1->indx;                /* set pO */
    for(i = 0; i < pO.size; i++){
        pO.data[i] = pE1->data[i];}
#if DISPLAYS
    printf("pO: ");
    ShowVector(&pO);
#endif
    for(i = 0; i < pO.size; i++){
        pO.data[i] = GFDiv(pE1->data[i], t);}
#if DISPLAYS
    printf("pO: ");
    ShowVector(&pO);
#endif

    if(Poly2Root(&vI, &pE)){
        printf("pE poly2root failed\n");
        return;
    }
#if DISPLAYS
    printf("vI: ");
    ShowVector(&vI);
#endif

    vO.size = vI.size;                  /* offsets */
    vX.size = vI.size;
    vY.size = vI.size;
    for(j = 0; j < vI.size; j++){
        i = vA2I.data[abInv[vI.data[j]]];
        vO.data[j] = i;
        vX.data[j] = vA.data[i];
        vY.data[j] = vU.data[i];
    }
#if DISPLAYS
    printf("vO: ");
    ShowVector(&vO);
    printf("vX: ");
    ShowVector(&vX);
    printf("vY: ");
    ShowVector(&vY);
#endif

    vV.size = vI.size;                  /* values */
    for(i = 0; i < vV.size; i++){
        t = GFMpy(GFSub(0, vX.data[i]), Poly2Val(&pO, abInv[vX.data[i]]));
        t = GFDiv(t, GFMpy(vY.data[i], Poly2Val(&pD, abInv[vX.data[i]])));
        vV.data[i] = t;
    }
    printf("vV: ");
    ShowVector(&vV);

    for(i = 0; i < vV.size; i++)        /* fix */
        vData.data[vO.data[i]] = GFSub(vData.data[vO.data[i]], vV.data[i]);
    printf("vData: ");
    ShowVector(&vData);

    vS.size = F;                        /* check syndromes */
    for(i = 0; i < F; i++)
        vS.data[i] = 0;
    for(j = 0; j < vData.size; j++){
        pS[j].size = F;
        for(i = 0; i < F; i++){
            pS[j].data[i] = GFMpy(pP[j].data[i], vData.data[j]);
            vS.data[F-1-i] = GFAdd(pS[j].data[i], vS.data[F-1-i]);
        }
    }
#if DISPLAYS
    printf("vS: ");
    ShowVector(&vS);
#endif
    t = 0;
    for(i = 0; i < vS.size; i++){
        if(vS.data[i] != 0){
            t = 1;
            break;
        }
    }
    if(t == 0){
        printf("no errors\n");
        return;
    }
    if(i != vS.size-1){
        printf("uncorrectable\n");
        for (i = 0; i < vV.size; i++)   /* unfix */
            vData.data[vO.data[i]] = GFAdd(vData.data[vO.data[i]], vV.data[i]);
        printf("vData: ");
        ShowVector(&vData);
        return;
    }
    t = vA2I.data[0];                   /* check for vA.data[] = 0 */
#if DISPLAYS                            /* fix single error at A2I[0] */
    printf("pP%x:", t);
    ShowVector(&pP[t]);
#endif
    vV.size = 1;
    vV.data[0] = GFDiv(vS.data[vS.size-1],pP[t].data[0]);
#if DISPLAYS
    printf("vV: ");
    ShowVector(&vV);
#endif
    vData.data[t] = GFSub(vData.data[t], vV.data[0]);
}
#endif

/*----------------------------------------------------------------------*/
/*      Fix     fix data                                                */
/*----------------------------------------------------------------------*/
static void Fix(void)
{
int i, j, k;
BYTE t;
    if((vE.size-1) > ((N-K)/2)){        /* check for e <= (N-K)/2 */
err0:
            printf("uncorrectable\n");
            return;
    }
    if(vP.size > K)                     /* check for deg(vP) >= K */
        goto err0;
    for(i = 0; i < vQ.size; i++){       /* check that rmdr = vQ == 0 */
        if(vQ.data[i] != 0)
            goto err0;
    }
    k = 0;                              /* find root indexes of vE */
    for(j = 0; j < N; j++){
        t = vE.data[0];                 /* calculate error indicator */
        for(i = 1; i < vE.size; i++)
            t = GFAdd(GFMpy(t, vA.data[j]), vE.data[i]);
        if(t != 0)                      /* if not error continue */
            continue;
        f[k++] = j;
    }
    if(k != vE.size-1)                  /* check for uncorrectable */
        goto err0;
    for(j = 0; j < vE.size; j++){
        t = vP.data[0];                 /* calculate value */
        for(i = 1; i < vP.size; i++)
            t = GFAdd(GFMpy(t, vA.data[f[j]]), vP.data[i]);
        vData.data[f[j]] = t;
    }
    printf("corrected %d errors\n", k);
}

/*----------------------------------------------------------------------*/
/*      PolyDiv(pPQuot, pPDvnd, pPDvsr) PQuot = PDvnd/PDvsr             */
/*                                      PDvnd = PDvnd%PDvsr             */
/*----------------------------------------------------------------------*/
static void PolyDiv(POLY *pPQuot, POLY *pPDvnd, POLY *pPDvsr)
{
int i, j;
BYTE t;

    pPQuot->size = pPDvnd->size - pPDvsr->size + 1;     /* poly divide */
    if(pPQuot->size <= 0){                              /*  handle quot == 0 */
        pPQuot->size = 1;
        pPQuot->data[0] = 0;
    }
    for(j = 0; j < pPQuot->size; j++){
        t = GFDiv(pPDvnd->data[j], pPDvsr->data[0]);
        pPQuot->data[j] = t;
        for(i = 0; i < pPDvsr->size; i++)
            pPDvnd->data[i+j] = GFSub(pPDvnd->data[i+j], GFMpy(t, pPDvsr->data[i]));
    }
}

/*----------------------------------------------------------------------*/
/*      MatrixMpy(pMDst, pMSrc0, pmSrc1) matrix multiply                */
/*----------------------------------------------------------------------*/
static void MatrixMpy(MATRIX *pMDst, MATRIX *pMSrc0, MATRIX *pMSrc1)
{
int i, j, k;
int drows, dcols, inner;
BYTE *pbd;                              /* ptr to dst */
BYTE *pbs0, *pbs1;                      /* ptr to src */

    drows = pMSrc0->nrows;              /* calc dest params */
    dcols = pMSrc1->ncols;
    inner = pMSrc0->ncols;              /* inner product size */

    pMDst->nrows = drows;               /* init dest */
    pMDst->ncols = dcols;
    memset(pMDst->data, 0, drows*dcols*sizeof(BYTE));

    pbd = pMDst->data;                  /* do the mpy */
    for(k = 0; k < drows; k++){
        for(j = 0; j < dcols; j++){
            pbs0 = &pMSrc0->data[k*inner];
            pbs1 = &pMSrc1->data[j];
            for(i = 0; i < inner; i++){
                *pbd = GFAdd(*pbd, GFMpy(*pbs0, *pbs1));
                pbs0 += 1;
                pbs1 += dcols;}
            pbd += 1;}}
}

/*----------------------------------------------------------------------*/
/*      MatrixInv(pMDst, pMSrc) invert matrix                           */
/*      assumes square matrix                                           */
/*----------------------------------------------------------------------*/
static int MatrixInv(MATRIX *pMDst, MATRIX *pMSrc)
{
int     i, j;
BYTE    *p0;                            /* generic pointers */
BYTE    *p1;
BYTE    *p2;
int     iNCol;                          /* # columns */
int     iNCol2;                         /* 2 * # columns */
int     iNCol21;                        /* 1+2*# columns */
BYTE    bMod;                           /* row modifier */
BYTE    bTemp;

    memset(abId, 0, sizeof(abId));      /* set up abId */
    abId[K] =  1;

    iNCol   = pMSrc->ncols;             /* set size related stuff */
    iNCol2  = iNCol+iNCol;
    iNCol21 = iNCol2+1;

    pMDst->nrows = iNCol;               /* set destination size */
    pMDst->ncols = iNCol2;              /*   for display */

/*      generate double-width augmented matrix */
/*      left side is copy of pMSrc, right side is Identity Matrix */

    p0  = pMSrc->data;
    p1  = pMDst->data;
    for(i = 0; i < iNCol; i++){
        memcpy(p1,  p0, iNCol*sizeof(BYTE)); /* copy a row of data */
        p0  += iNCol;
        p1  += iNCol;
        memcpy(p1,  &abId[K-i], iNCol*sizeof(BYTE)); /* add ID  part */
        p1  += iNCol;}

/*      normalize according to left side */
/*      results in inverse matrix in right size */

    #if DISPLAYI
        printf("start\n");
        ShowMatrix(pMDst);
    #endif

    for(j = 0; j < iNCol; j++){         /* working at [j][j] */

/*      find 1st non-zero in current column */

        p0  = pMDst->data+j*iNCol21;    /* p0 = starting ptr */
        p1  = p0;                       /* p1 = scanning ptr */
        p2  = pMDst->data+iNCol2*iNCol; /* p2 = end of  matrix */
        while(0 == *p1){
            p1  += iNCol2;
            if(p1 >= p2){               /* return if bad matrix */
                return(1);}}

        bMod = *p1;                     /* set divisor for row */

/*      swap rows if needed */

        if(p0 != p1){
            p0  -= j;                   /* backup to beg of rows */
            p1  -= j;
            for(p2  = p0+iNCol2; p0 != p2; p0++, p1++){
                bTemp = *p0;
                *p0 = *p1;
                *p1 = bTemp;}
            #if DISPLAYI
                printf("swapped rows\n");
                ShowMatrix(pMDst);
            #endif
            ;}


/*      divide row to produce a one */

        p0  = pMDst->data+j*iNCol2;     /* p0 = ptr to  start of row */
        for(p2  = p0+iNCol2; p0 != p2; p0++){
            *p0 = GFDiv(*p0, bMod);}

        #if DISPLAYI
            printf("divided row\n");
            ShowMatrix(pMDst);
        #endif

/*      subtract multiple of this row from other rows */
/*      to create a column of zeroes */
/*      (except for this row,column) */

        for(i = 0; i < iNCol; i++){
            if(i == j)                  /* skip if current row */
                continue;
            p0  = pMDst->data+j*iNCol2; /* p0 = ptr to  current row */
            p1  = pMDst->data+i*iNCol2; /* p1 = ptr to  target row */
            bMod = *(p1+j);             /* bMod = current row mpyr */
            for(p2  = p0+iNCol2; p0 != p2; p0++, p1++)
                *p1 = GFSub(*p1, GFMpy(*p0, bMod));}
#if DISPLAYI
        printf("zeroed columns\n");
        ShowMatrix(pMDst);
#endif
        ;}

/*      now copy right side of matrix to left side */

    p0  = pMDst->data;                  /* p0 = Dst ptr */
    p1  = p0+iNCol;
    for(j = 0; j < iNCol; j++){
        p2  = p0+iNCol;
        while(p0 != p2){
            *p0++ = *p1++;}
        p1  += iNCol;}

    pMDst->ncols = iNCol;               /* set proper ncols */
#if DISPLAYI
    printf("shifted left\n");
    ShowMatrix(pMDst);
#endif
    return(0);
}

/*----------------------------------------------------------------------*/
/*      Poly2Root(pVDst, pPSrc)         find roots of poly from vA      */
/*----------------------------------------------------------------------*/
static int Poly2Root(VECTOR *pVDst, VECTOR *pPSrc)
{
int     i,j;
int     iA;                             /* index to vA */
int     iDst;                           /* index to pVDst */
BYTE    bRoot;                          /* current root */
BYTE    bSum;                           /* current sum */

    pVDst->size = pPSrc->size-1;        /* set dest size */

    if(!pVDst->size)                    /* exit if null */
        return(0);

    iA = 0;
    iDst   = 0;
    for(j = 0; j < N;  j++){
        bRoot = vA.data[iA++];
        if(bRoot == 0)
            continue;
        bRoot = GFDiv(1, bRoot);
        bSum = 0;                       /* sum up terms */
        for(i = 0; i < pPSrc->size; i++){
            bSum = GFMpy(bSum, bRoot);
            bSum = GFAdd(bSum, pPSrc->data[i]);}

        if(!bSum){                      /* if a root */
            if(iDst == pVDst->size){    /*    exit if too many roots */
                return(1);}
            pVDst->data[iDst] = bRoot;  /*    add root */
            iDst++;
        }
    }

    if(iDst != pVDst->size)             /* exit if not enough roots */
        return(1);

    return(0);
}

/*----------------------------------------------------------------------*/
/*      Root2Poly(pPDst, pVSrc)         convert roots into polynomial   */
/*----------------------------------------------------------------------*/
static void Root2Poly(VECTOR *pPDst, VECTOR *pVSrc)
{
int i, j;

    pPDst->size = pVSrc->size+1;
    pPDst->data[0] = 1;
    memset(&pPDst->data[1], 0, pVSrc->size*sizeof(BYTE));
    for(j = 0; j < pVSrc->size; j++){
        for(i = j; i >= 0; i--){
            pPDst->data[i+1] = GFSub(pPDst->data[i+1],
                    GFMpy(pPDst->data[i], pVSrc->data[j]));}}
}

/*----------------------------------------------------------------------*/
/*      Poly2Val(pPoly, val)        evaluates polynomial(val)           */
/*----------------------------------------------------------------------*/
static BYTE Poly2Val(VECTOR *pPoly, BYTE val)
{
int s = pPoly->size-1;
int i;
BYTE v = 0;
    for(i = 0; i < pPoly->size; i++)
        v = GFAdd(v, GFMpy(GFPow(val, i), pPoly->data[s-i]));
    return v;
}

/*----------------------------------------------------------------------*/
/*      GFAdd(b0, b1)           b0+b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFAdd(BYTE b0, BYTE b1)
{
    return(b0^b1);
}

/*----------------------------------------------------------------------*/
/*      GFSub(b0, b1)           b0-b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFSub(BYTE b0, BYTE b1)
{
    return(b0^b1);
}

/*----------------------------------------------------------------------*/
/*      GFMpy(b0, b1)           b0*b1       using logs                  */
/*----------------------------------------------------------------------*/
static BYTE GFMpy(BYTE byt0, BYTE byt1)
{
    if(byt0 == 0 || byt1 == 0)
        return(0);
    return(abExp[(int)abLog[byt0]+(int)abLog[byt1]]);
}

/*----------------------------------------------------------------------*/
/*      GFDiv(b0, b1)           b0/b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFDiv(BYTE b0, BYTE b1)
{
    if(b1 == 0){
        printf("divide by zero\n");
        return(0);
    }
    if(b0 == 0)
        return(0);
    return(abExp[15+(int)abLog[b0]-(int)abLog[b1]]);
}

/*----------------------------------------------------------------------*/
/*      GFPow(b0, b1)           b0^b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFPow(BYTE b0, int i1)
{
    if(i1 == 0)
        return (1);
    if(i1 < 0){
        printf("negative pow\n");
        i1 += 15;
    }
    if(b0 == 0)
        return (0);
    return(abExp[(abLog[b0]*i1)%(15)]);
}

/*----------------------------------------------------------------------*/
/*      GFMpy0(b0,b1)           b0*b1       using low level math        */
/*----------------------------------------------------------------------*/
static BYTE GFMpy0(BYTE b0, BYTE b1)
{
int i;
int product;
    product = 0;
    for(i = 0; i < 4; i++){
        product <<= 1;
        if(product & 0x10u){
            product ^= iGF;}
        if(b0 & 0x8u){
            product ^= b1;}
        b0 <<= 1;}
    return((BYTE)product);
}

/*----------------------------------------------------------------------*/
/*      InitGF  Initialize Galios Stuff                                 */
/*----------------------------------------------------------------------*/
static void InitGF(void)
{
int i;
#if SYN
int j, k;
#endif
BYTE t;

    t = 1;                              /* init abExp[] */
    for(i = 0; i < 32; i++){
        abExp[i] = t;
        t = GFMpy0(t, bAlpha);
    }

    abLog[0] = 0xffu;                   /* init abLog[] */
    for(i = 0; i < 15; i++)
        abLog[abExp[i]] = i;

    abInv[0] = 0xffu;                   /* init abInv[] */
    for(i = 1; i < 16; i++)
        abInv[i] = GFDiv(1,i);

    abNeg[0] = 0;                       /* init abNeg[] */
    for(i = 1; i < 16; i++)
        abNeg[i] = 16-i;

    vA.size = N;                        /* init vA */
#if 1
    memcpy(vA.data, "\x7\xf\x0\x1\x5\xd\x9\x3\xa\x2\x8\xb\x4\xc\x6\xe", 16);
#else
    memcpy(vA.data, "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xa\xb\xc\xd\xe\xf", 16);
#endif
#if DISPLAYX
    printf("vA: ");
    ShowVector(&vA);
#endif

    Vandermonde();                      /* generate Vandermonde matrix */

    LagrangeI(&vLagrK, K);              /* generate Lagrange */
    LagrangeI(&vLagrN, N);

#if SYN
    vA2I.size = vA.size;
    vA2I.data[0] = vA.size;
    for(i = 0; i < vA2I.size; i++)
        vA2I.data[vA.data[i]] = i;
    vB.size = vA.size-1;
    for(j = 0; j < vA.size; j++){
        k = 0;
        for(i = 0; i < vA.size; i++){
            if(i == j){
                continue;
            }
            vB.data[k++] = vA.data[i];
        }
        Root2Poly(&pL[j], &vB);
    }
    vL.size = vA.size;
    for(i = 0; i < vA.size; i++)
        vL.data[i] = Poly2Val(&pL[i], vA.data[i]);
#if DISPLAYX
    printf("vL: ");
    ShowVector(&vL);
#endif
    vU.size = vA.size;
    for(i = 0; i < vA.size; i++)
        vU.data[i] = GFDiv(1, vL.data[i]);
#if DISPLAYX
    printf("vU: ");
    ShowVector(&vU);
#endif
    for(j = 0; j < vA.size; j++){
        pP[j].size = F;
        for(i = 0; i < F; i++){
            pP[j].data[i] = GFPow(vA.data[j], i);
        }
        for(i = 0; i < F; i++){
            pP[j].data[i] = GFMpy(pP[j].data[i], vU.data[j]);
        }
#if DISPLAYX
        printf("pP%1x ", j);
        ShowVector(&pP[j]);
#endif
    }
#endif
}

/*----------------------------------------------------------------------*/
/*      ShowVector                                                      */
/*----------------------------------------------------------------------*/
static void ShowVector(VECTOR *pVSrc)
{
int i;
    for(i = 0; i < pVSrc->size; i++)
        printf(" %1x", pVSrc->data[i]);
    printf("\n");
}

/*----------------------------------------------------------------------*/
/*      ShowEuclid
/*----------------------------------------------------------------------*/
static void ShowEuclid(EUCLID *pESrc)
{
int i;
    for(i = 0; i < pESrc->indx; i++)
        printf(" %1x", pESrc->data[i]);
    printf("|");
    for( ; i < pESrc->size; i++)
        printf("%1x ", pESrc->data[i]);
    printf("\n");
}

/*----------------------------------------------------------------------*/
/*      ShowMatrix                                                      */
/*----------------------------------------------------------------------*/
static void ShowMatrix(MATRIX *pMSrc)
{
int iM;
int i, j;
    iM = 0;
    for(i = 0; i < pMSrc->nrows; i++){
        for(j = 0; j < pMSrc->ncols; j++)
            printf(" %1x", pMSrc->data[iM++]);
        printf("\n");
    }
}

/*----------------------------------------------------------------------*/
/*      Conrs get console response - emulates _cgets()                  */
/*----------------------------------------------------------------------*/
static int Conrs(void)
{
int i;
    memset(abUser, 0, sizeof(abUser));  /* get a line */
    abUser[0] = sizeof(abUser)-2;
    fgets(abUser+2, sizeof(abUser)-2, stdin);
    abUser[1] = (int)(strlen(&abUser[2])-1);
    i = abUser[1];
    abUser[2+i] = 0;
    return(i);
}

/*----------------------------------------------------------------------*/
/*      main                                                            */
/*----------------------------------------------------------------------*/
int main()
{
    DoUser();
    return(0);
}

/*----------------------------------------------------------------------*/
/*      DoUser  do user stuff                                           */
/*----------------------------------------------------------------------*/
static void DoUser(void)
{
int i, j, u0, u1;

    printf("eccr16z 1.0\n");

    while(1){
        printf("Galios Field Generators and Alpha's\n");
        printf("   ");
        for(i = 0; i < 3; i++){
            printf("%7d", i);}
        printf("\n");
        u0 = 0;
        for(j = 0; j < 1; j++){
            printf("%02d:", u0/2);
            for(i = 0; i < 3; i++){
                printf(" %03x  %01x", aiGA[u0], aiGA[u0+1]);
                u0 += 2;}
            printf("\n");}
        printf("Select Galios Field Generator (0-2):\n");
        if(!Conrs())return;
        sscanf(&abUser[2], "%d", &u0);
        if(u0 > 2)continue;
        u0  <<= 1;
        iGF    = aiGA[u0];
        bAlpha  = (BYTE)(aiGA[u0+1]);
        printf(    "    %03x  %01x\n", iGF, bAlpha);
        break;}
    while(1){
        printf("Enter # of parity values (1->%d):\n", MAXPAR);
        if(!Conrs())return;
        sscanf(&abUser[2], "%d", &u0);
        if(u0 < 1 || u0 > MAXPAR)continue;
        break;
    }
    while(1){
        printf("Enter # data values (1->%d):\n", MAXDAT-u0);
        if(!Conrs())return;
        sscanf(&abUser[2], "%d", &u1);
        if(u1 < 1 || u1 > (MAXDAT-u0))continue;
        K = u1;
        N = u0+u1;
        F = N-K;
        E = F/2;
        vData.size = N;
        break;
    }
    InitGF();
DoUser0:
    printf("vData:\n");
    ShowVector(&vData);
#if 0
    printf("vPoly:\n");
    ShowVector(&vPoly);
#endif

/*      get user command */
#if FOURIER
    printf("Enter Change Zero Trans Utrns Ncode Scode Vcode Xvcode Fix Quit: ");
#else
    printf("Enter Change Zero Ncode Scode Vcode Xvcode Fix Quit: ");
#endif
    if(0 == Conrs())
        goto DoUser0;
    printf("\n");
    switch(abUser[2] & 0x5f){
      case 'E':                         /* enter data */
        for(i = 0; i < vData.size; i++){
            printf("%2x = ", i);
            if(!Conrs())break;
            sscanf(&abUser[2], "%x", &u0);
            vData.data[i] = u0;
        }
        break;
      case 'C':                         /* change data */
        while(1){
            printf("Enter offset (0-%x): ", vData.size-1);
            if(!Conrs())break;
            sscanf(&abUser[2], "%x", &u0);
            if(u0 >= vData.size)continue;
            printf("Enter value: (0->f): ");
            if(!Conrs())break;
            sscanf(&abUser[2], "%x", &u1);
            if(u1 >= 0x10)continue;
            vData.data[u0] = u1;
        }
        break;
#if 0
      case 'P':                         /* enter poly */
        vPoly.size = K;
        for(i = 0; i < K; i++){
            printf(" x^%1d c: ", K-1-i);
            if(!Conrs())break;
            sscanf(&abUser[2], "%d", &u0);
            vPoly.data[i] = u0;
        }
        EncodeP();
        break;
#endif
      case 'Z':                         /* zero data array */
        for(i = 0; i < N; i++)
            vData.data[i] = 0;
        Encode();
        break;
#if FOURIER
      case 'T':
        Fourier();
        break;
      case 'U':
        FourierI();
        break;
#endif
      case 'N':                         /* encode */
        Encode();
        break;
      case 'S':                         /* encode systematic */
        EncodeS();
        break;
      case 'V':                         /* encode vandermonde */
        EncodeV();
        break;
      case 'X':                         /* undo encode vandermonde */
        UnEncodeV();
        Encode();
        break;
      case 'F':                         /* fix (DecodeBW) */
#if BW
        DecodeBW();
        Fix();
#endif
#if GAO
        DecodeG();
        Fix();
#endif
#if SYN
        DecodeS();
#endif
        break;
      case 'Q':                         /* quit */
        return;
      default:
        break;
    }
    goto DoUser0;
}
