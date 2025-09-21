/*----------------------------------------------------------------------*/
/*      bch1022.c       bch(1022,992)                                   */
/*                                                                      */
/*      Jeff Reid       2021JAN09 11:30                                 */
/*----------------------------------------------------------------------*/
#include <intrin.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned char      BYTE;
typedef unsigned short     WORD;
typedef unsigned long      DWORD;
typedef unsigned long long QWORD;
typedef          short     SWORD;

/* GF(2^10) = x^10 + x^3 + 1 */

#define POLY 0x409

/* GF(2^10) primitive */
#define ALPHA 0x0002

/* BCH(1022, 992), BN, BK are muliples of 8 */
#define BN 1024
#define BK 992

/* number of suyndromes */
#define NSYN 6

/*                              ** if != 0, use extended euclid algorithm */
#define EEUCLID 0

/*                              ** display euclid stuff */
#define DISPLAYE 0

/*                              ** display stuff */
#define DISPLAY 0

typedef struct{                 /* vector structure */
    WORD  size;
    WORD  data[31];
}VECTOR;

#if EEUCLID
typedef struct{                 /* euclid structure */
    WORD  size;                 /* # of data words */
    WORD  indx;                 /* index to right side */
    WORD  data[NSYN+2];         /* left and right side data */
}EUCLID;
#endif

/*                              ** if QP != 0, use queryperformance for timer */
#define QP 1

#if QP
#include <math.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
typedef LARGE_INTEGER LI64;
#else
#include <time.h>
#endif

/*----------------------------------------------------------------------*/
/*      data                                                            */
/*----------------------------------------------------------------------*/
static WORD gflog2[1024];       /* gflog2 table */
static WORD gfexp2[1024];       /* gfexp2 table */

static WORD minply[3];          /* minimum polynomials */
static WORD mincnt;             /* # of minimum polymials */
static BYTE minplyf[1024];      /* minimum polynomial flags */

static DWORD polytbl[256];      /* encode poly table */

static WORD syntbl[256][NSYN];  /* syndrome table */

static BYTE msg[BN/8];          /* encoded message, last 2 bits unused */

#if EEUCLID
static EUCLID   E0;             /* used by GenpErrors (extended Euclid) */
static EUCLID   E1;
#else
static VECTOR   vB;             /* used by GenpErrors (Berleykamp Massey) */
static VECTOR   vC;
static VECTOR   vT;
static VECTOR   vBx;
#endif

static VECTOR   vSyndromes;
static VECTOR   pErrors;
static VECTOR   pLambda;
static VECTOR   vLocators;
static VECTOR   vOffsets;

#if QP
static LI64     liQPFrequency;  /* cpu counter values */
static LI64     liStartTime;
static LI64     liStopTime;
static double   dQPFrequency;
static double   dStartTime;
static double   dStopTime;
static double   dElapsedTime;
#else
static clock_t ctTimeStart;     /* clock values */
static clock_t ctTimeStop;
#endif

/*----------------------------------------------------------------------*/
/*      code                                                            */
/*----------------------------------------------------------------------*/
static void GenMinPoly(void);
static void GenPolyTbl(void);
static void GenSynTbl(void);
static void Encode(void);
static void GenSyndromes(void);
static void GenpErrors(void);
static void GenOffsets(void);
static void FixErrors(void);
static int  Poly2Root(VECTOR *, VECTOR *);
static void GFInit(void);
static WORD GFPwr(WORD, WORD);
static WORD GFMpy(WORD, WORD);
static WORD GFDiv(WORD, WORD);
#define GFAdd(a, b) (a^b)
#define GFSub(a, b) (a^b)
static void ShowVector(VECTOR *);
#if EEUCLID
static void ShowEuclid(EUCLID *);
#endif

/*----------------------------------------------------------------------*/
/*      main                                                            */
/*----------------------------------------------------------------------*/
main()
{
WORD i;
    GFInit();                   /* init GF(2^10) tables */
    GenMinPoly();               /* generate minimum polynomials */
    GenPolyTbl();               /* generate 192 bit poly table */
    GenSynTbl();                /* generate syndrome table */
    for(i = 0; i < (BK/8); i++) /* generate test message */
        msg[i] = (BYTE)i;
#if QP
    QueryPerformanceFrequency(&liQPFrequency);
    dQPFrequency = (double)liQPFrequency.QuadPart;
    QueryPerformanceCounter(&liStartTime);
#else
    ctTimeStart = clock();
#endif
    Encode();
#if QP
    QueryPerformanceCounter(&liStopTime);
    dStartTime = (double)liStartTime.QuadPart;
    dStopTime  = (double)liStopTime.QuadPart;
    dElapsedTime = (dStopTime - dStartTime) / dQPFrequency;
    printf("# of seconds %f\n", dElapsedTime);
#else
    ctTimeStop = clock();
    printf("# of ticks %u\n", ctTimeStop - ctTimeStart);
#endif
      msg[  0] ^= 0x80;         /* test 3 error bit case */
      msg[ 64] ^= 0x10;
      msg[127] ^= 0x08;
#if QP
    QueryPerformanceCounter(&liStartTime);
#else
    ctTimeStart = clock();
#endif
    GenSyndromes();
    GenpErrors();
    GenOffsets();
    FixErrors();
#if QP
    QueryPerformanceCounter(&liStopTime);
    dStartTime = (double)liStartTime.QuadPart;
    dStopTime  = (double)liStopTime.QuadPart;
    dElapsedTime = (dStopTime - dStartTime) / dQPFrequency;
    printf("# of seconds %f\n", dElapsedTime);
#else
    ctTimeStop = clock();
    printf("# of ticks %u\n", ctTimeStop - ctTimeStart);
#endif

#if DISPLAY
    printf("pErrors:  ");
    ShowVector(&pErrors);
    printf("vLocators:");
    ShowVector(&vLocators);
    printf("vOffsets: ");
    ShowVector(&vOffsets);
#endif
    GenSyndromes();
    for(i = 0; i < NSYN; i++)
        if(vSyndromes.data[i] != 0)
            break;
    if(i == NSYN)
        printf("passed\n");
    else
        printf("failed\n");
    return 0;
}

/*----------------------------------------------------------------------*/
/*      GenMinPoly      generate 3 mininum polynomials                  */
/*----------------------------------------------------------------------*/
static void GenMinPoly(void)
{
WORD poly;                      /* test polynomial */
WORD sum;                       /* sum, looking for zeroes */
WORD apwr;                      /* alpha to power */
WORD i,j;

    /* find 3 minimum polynomials for 5 powers of 2 */
    for(i = 1; i < NSYN; i++){
        apwr = GFPwr(2,i);
        for(poly = 0x401; poly <= 0x7ff ; poly++){
            sum = 0;
            for(j = 0; j <= 10; j++){
                if(poly&(1<<j))
                    sum ^= GFPwr(apwr,j);
            }
            if(sum == 0){
                if(minplyf[poly-0x400] != 0)
                    continue;
                minplyf[poly-0x400] = 1;
                minply[mincnt++] = poly;
                break;
            }
        }
    }
}

/*----------------------------------------------------------------------*/
/*      GenPolyTbl     generate 32 bit polynomial table                 */
/*----------------------------------------------------------------------*/
static void GenPolyTbl(void)
{
__m128i p;                      /* sub-products and product */
DWORD poly;                     /* 32 bit polynomial */
DWORD t;                        /* for table generation */
DWORD b, q, i;                  /* byte, quotient bit, i */

/*  generate poly */
    p.m128i_u64[0] = minply[0];
    p.m128i_u64[1] = minply[1];
    p              = _mm_clmulepi64_si128(p, p, 0x01);
    p.m128i_u64[1] = minply[2];
    p = _mm_clmulepi64_si128(p, p, 0x01);
    poly = p.m128i_u32[0]<<2;
/* generate polytbl */
    for(b = 0x000; b < 0x100; b++){
        t = b<<24;
        for(i = 0; i < 8; i++){
            q = t>>31;
            t = t<<1;
            if(q != 0){
                t ^= poly;
            }
        }
        polytbl[b] = t;
    }
}

/*----------------------------------------------------------------------*/
/*      GenSynTbl   generate syndrome table                             */
/*----------------------------------------------------------------------*/
static void GenSynTbl(void)
{
WORD i, j, k, apwr, sum;
    for(k = 0; k < 0x100; k++){
        for(j = 0; j < NSYN; j++){
            apwr = GFPwr(2, j+1);
            sum = 0;
            for(i = 0; i < 8; i++){
                if(k & (1<<i))
                    sum ^= GFPwr(apwr, i);
            }
            syntbl[k][j] = sum;
        }
    }
}

/*----------------------------------------------------------------------*/
/*      Encode                                                          */
/*----------------------------------------------------------------------*/
static void Encode(void)
{
DWORD r;                        /* remainder */
QWORD i, j;
    r = 0;
    for(j = 0; j < BK/8; j += 1){
        i = (r>>24)^msg[j];     /* table index */
        r = r<<8;
        r ^= polytbl[i];
    }
    msg[(BK/8)+ 0] = (BYTE)(r>>24);
    msg[(BK/8)+ 1] = (BYTE)(r>>16);
    msg[(BK/8)+ 2] = (BYTE)(r>> 8);
    msg[(BK/8)+ 3] = (BYTE)(r>> 0);
}

/*----------------------------------------------------------------------*/
/*      GenSyndromes                                                    */
/*----------------------------------------------------------------------*/
static void GenSyndromes(void)
{
WORD i, j, apwr;
    vSyndromes.size = NSYN;
    memset(vSyndromes.data, 0, NSYN*sizeof(WORD));
    for(i = 0; i < NSYN; i++){
        /* do byte at a time until last byte */
        apwr = GFPwr(2, (i+1)<<3);
        for(j = 0; j < (BN/8)-1; j++){
            vSyndromes.data[i] = GFMpy(vSyndromes.data[i], apwr) ^
                                 syntbl[msg[j]][i];
        }
        /* do last 6 bits */
        apwr = GFPwr(2, i+1);
        for(j = BN-8; j < BN-2; j++){
            vSyndromes.data[i] = GFMpy(vSyndromes.data[i], apwr) ^
                                 ((msg[j/8]>>(7-(j%8)))&1);
        }
    }
}

#if EEUCLID
/*----------------------------------------------------------------------*/
/*      GenpErrors     generate pErrors via Euclid division algorithm   */
/*----------------------------------------------------------------------*/
static void GenpErrors(void)
{
/* R[] is msb first | A[] is msb last (reversed) */
EUCLID *pED;                            /* R[i-2] | A[i-1] */
EUCLID *pER;                            /* R[i-1] | A[i-2] */
EUCLID *pET;                            /* temp */
int     i, j;
WORD    bME;                            /* max errors possible */
WORD    bQuot;                          /* quotient */

/*      E0 = initial ED: E0.R[-1] = x^MAXERR, E0.A[0] = 1 */
    E0.size = vSyndromes.size+2;
    E0.indx = vSyndromes.size+1;
    E0.data[0] = 1;
    memset(&E0.data[1], 0, vSyndromes.size*sizeof(WORD));
    E0.data[E0.indx] = 1;
    pED = &E0;

/*      E1 = initial ER: E1.R[0] = syndrome polynomial, E1.A[-1] = 0 */
    E1.size = vSyndromes.size+2;
    E1.indx = vSyndromes.size+1;
    E1.data[0] = 0;
    for(i = 1; i < E1.indx; i++){
        E1.data[i] = vSyndromes.data[vSyndromes.size-i];}
    E1.data[E1.indx] = 0;
    pER = &E1;

/*      init bME */
    bME = vSyndromes.size/2;

/*      Euclid algorithm */

    while(1){                           /* while degree ER.R > max errors */
#if DISPLAYE
        printf("ED: ");
        ShowEuclid(pED);
        printf("ER: ");
        ShowEuclid(pER);
#endif
        while((pER->data[0] == 0) &&    /* shift dvsr left until msb!=0 */
              (pER->indx != 0)){        /*  or fully shifted left */
            pER->indx--;
            memcpy(&pER->data[0], &pER->data[1], (pER->size-1)*sizeof(WORD));
            pER->data[pER->size-1] = 0;}

        if(pER->indx <= bME){           /* if degree ER.R[] <= bME, break */
            break;}

        while(1){                       /* while more sub-steps */
            if(pED->data[0]){           /*   if ED.R[] msb!=0, update ED, ER */
                bQuot = GFDiv(pED->data[0], pER->data[0]); /* Q=ED.R[msb]/ER.R[msb] */
                for(i = 0; i < pER->indx; i++){            /* ED.R[]=ED.R[]-Q*ER.R[] */
                    pED->data[i] = GFSub(pED->data[i], GFMpy(bQuot, pER->data[i]));}
                for(i = pED->indx; i < pER->size; i++){    /* ER.A[]=ER.A[]-Q*ED.A[] */
                    pER->data[i] = GFSub(pER->data[i], GFMpy(bQuot, pED->data[i]));}}
            if(pED->indx == pER->indx){ /*   if sub-steps done, break */
                break;}
            pED->indx--;                /*   shift ED left */
            memcpy(&pED->data[0], &pED->data[1], (pED->size-1)*sizeof(WORD));
            pED->data[pED->size-1] = 0;}

        pET = pER;                      /* swap ED, ER */
        pER = pED;
        pED = pET;}

    pErrors.size = pED->size-pED->indx; /* set pErrors.size */

    if((pER->indx) >= pErrors.size){    /*  if degree ER.R too high */
        printf("GenpErrors remainder.size >= errors.size\n");
        goto fail0;}

#if 0
    j = pErrors.size - 1;       /* right shift ER if Omega has leading zeroes */
    while(pER->indx < j){
        pER->indx++;
        for(i = pER->size-1; i;){
            i--;
            pER->data[i+1] = pER->data[i];}
        pER->data[0] = 0;}
#if DISPLAYE
    printf("EX: ");
    ShowEuclid(pER);
#endif
#endif

/*      pErrors = ED.A[] without unreversing = Lambda reversed */
    j = pED->indx;
    for(i = 0; i < pErrors.size; i++){
        pErrors.data[i] = pED->data[j];
        j++;}

#if DISPLAYE
    printf("pErrors (e):    ");
    ShowVector(&pErrors);
#endif

/*      Make most significant coef pErrors == 1 (divide by it) */
    bQuot = pErrors.data[0];
    if(bQuot == 0){
        printf("GenpErrors most sig coef of pErrors == 0\n");
        pLambda.size = 1;
        pLambda.data[0] = 1;
        goto fail0;}
    for(i = 0; i < pErrors.size; i++){
        pErrors.data[i] = GFDiv(pErrors.data[i], bQuot);}
#if DISPLAYE
    printf("pErrors (E):    ");
    ShowVector(&pErrors);
#endif

/*      Find roots of pErrors (if fail, then set for no roots) */

    if(Poly2Root(&vLocators, &pErrors)){    /* solve error poly */
        printf("GenpErrors poly2root failed \n");
fail0:
        pErrors.size = 1;                   /* handle no root case */
        pErrors.data[0] = 1;
        vLocators.size = 0;}
}

#else

/*----------------------------------------------------------------------*/
/*      GenpErrors      generate pErrors via Berklekamp Massey          */
/*      note poly most signifcant index == 0                            */
/*----------------------------------------------------------------------*/
static void GenpErrors(void)
{
WORD i, j, n;
WORD L, m;
WORD b, d;
WORD db;

    b = 1;                              /* discrepancy when L last updated */
    L = 0;                              /* number of errors */
    m = 1;                              /* # iterations since L last updated */
    vB.size    = 1;
    vB.data[0] = 1;
    vC.size    = 1;
    vC.data[0] = 1;

    for(n = 0; n < vSyndromes.size; n++){
        if(n&1){                        /* BCH only, if odd step, d == 0 */
            m += 1;
            continue;}
        d = vSyndromes.data[n];         /* calculate discrepancy */
        for(i = 1; i <= L; i++){
            d = GFAdd(d, GFMpy(vC.data[(vC.size - 1)- i], vSyndromes.data[n-i]));}
        if(d == 0){                     /* if 0 increment m, continue */
            m += 1;
            continue;}
        vT.size = vC.size;              /* vT = vC */
        memcpy(vT.data, vC.data, vC.size*sizeof(WORD));
        db = GFDiv(d,b);                /* db = (d/b) */
        vBx.size = vB.size+m;           /* Bx = x^m B */
        memcpy(vBx.data, vB.data, vB.size*sizeof(WORD));
        memset(&vBx.data[vB.size], 0, m*sizeof(WORD));
        for(i = 0; i < vBx.size; i++){  /* Bx *= db */
            vBx.data[i] = GFMpy(vBx.data[i], db);}
        j = vBx.size - vC.size;         /* right shift vBx or vC */
        if(((SWORD)j) > 0){
            for(i = vBx.size; i > j; ){
                i--;
                vC.data[i] = vC.data[i-j];}
            memset(vC.data, 0, j*sizeof(WORD));
            vC.size += j;}
        else if(((SWORD)j) < 0){
            j = -j;
            for(i = vC.size; i > j; ){
                i--;
                vBx.data[i] = vBx.data[i-j];}
            memset(vBx.data, 0, j*sizeof(WORD));
            vBx.size += j;}
        for(i = 0; i < vC.size; i++){   /* C -= Bx */
            vC.data[i] = GFSub(vC.data[i], vBx.data[i]);}
        if(n < 2*L){                    /* if L not increasing */
            m += 1;
            continue;}
        vB.size = vT.size;              /*   B = T */
        memcpy(vB.data, vT.data, vT.size*sizeof(WORD));
        L = n + 1 - L;                  /* update L */
        b = d;                          /* update b */
        m = 1;}                         /* reset m */

    pErrors.size = vC.size;             /* pErrors = reversed VC */
    for(i = 0; i < vC.size; i++)
        pErrors.data[i] = vC.data[vC.size-1-i];

    if(Poly2Root(&vLocators, &pErrors)){    /* solve error poly */
        printf("GenpErrors poly2root failed \n");
        pErrors.size = 1;                   /* handle no root case */
        pErrors.data[0] = 1;
        vLocators.size = 0;}
}

#endif

/*----------------------------------------------------------------------*/
/*      GenOffsets                                                      */
/*----------------------------------------------------------------------*/
static void GenOffsets(void)
{
WORD i;
    vOffsets.size = vLocators.size;
    for(i = 0; i < vLocators.size; i++){
        vOffsets.data[i] = BN-3-gflog2[vLocators.data[i]];
    }
}

/*----------------------------------------------------------------------*/
/*      FixErrors                                                       */
/*----------------------------------------------------------------------*/
static void FixErrors()
{
WORD i;
    for(i = 0; i < vOffsets.size; i++)
        msg[vOffsets.data[i]/8] ^= (WORD)(0x80u>>(vOffsets.data[i]%8));
}

/*----------------------------------------------------------------------*/
/*      Poly2Root(pVDst, pPSrc)         find roots of poly              */
/*----------------------------------------------------------------------*/
static int Poly2Root(VECTOR *pVDst, VECTOR *pPSrc)
{
WORD    bLcr;                           /* current locator */
WORD    bSum;                           /* current sum */
WORD    bV;                             /* index to pVDst */
WORD    i,j;

    pVDst->size = pPSrc->size-1;        /* set dest size */

    if(!pVDst->size)                    /* exit if null */
        return(0);

    bV   = 0;
    bLcr = 1;
    for(j = 0; j < BN;  j++){
        bSum = 0;                       /* sum up terms */
        for(i = 0; i < pPSrc->size; i++){
            bSum = GFMpy(bSum, bLcr);
            bSum = GFAdd(bSum, pPSrc->data[i]);}

        if(!bSum){                      /* if a root */
            if(bV == pVDst->size){      /*    exit if too many roots */
                return(1);}
            pVDst->data[bV] = bLcr;     /*    add locator */
            bV++;}

        bLcr = GFMpy(bLcr, ALPHA);}     /* set up next higher alpha */

    if(bV != pVDst->size)               /* exit if not enough roots */
        return(1);

    return(0);
}

/*----------------------------------------------------------------------*/
/*      GFInit          init GF(2^10) tables                            */
/*----------------------------------------------------------------------*/
static void GFInit(void)
{
WORD *p0, *p1;
WORD d0;
    d0 = 1;                             /* init gfexp2 table */
    p0 = gfexp2;
    for(p1 = p0+1024; p0 < p1;){
        *p0++ = d0;
        d0 <<= 1;
        if(d0 & 0x400)
            d0 ^= POLY;}

    p0 = gfexp2;                        /* init gflog2 table */
    p1 = gflog2;
    *p1 = 0;
    for(d0 = 0; d0 < 1023; d0 += 1)
        *(p1+*p0++) = d0;
}

/*----------------------------------------------------------------------*/
/*      GFPwr - exponentiate                                            */
/*----------------------------------------------------------------------*/
static WORD GFPwr(WORD m0, WORD m1)
{
    return gfexp2[(WORD)((gflog2[m0]*(DWORD)m1)%1023)];
}

/*----------------------------------------------------------------------*/
/*      GFMpy                                                           */
/*----------------------------------------------------------------------*/
static WORD GFMpy(WORD m0, WORD m1)   /* multiply */
{
int m2;
    if(0 == m0 || 0 == m1)
        return(0);
    m2 = gflog2[m0] + gflog2[m1];
    if(m2 > 1023)
        m2 -= 1023;
    return(gfexp2[m2]);
}

/*----------------------------------------------------------------------*/
/*      GFDiv                                                           */
/*----------------------------------------------------------------------*/
static WORD GFDiv(WORD m0, WORD m1)   /* divide */
{
int m2;
    if(0 == m0)
        return(0);
    m2 = gflog2[m0] - gflog2[m1];
    if(m2 < 0)
        m2 += 1023;
    return(gfexp2[m2]);
}

/*----------------------------------------------------------------------*/
/*      ShowVector                                                      */
/*----------------------------------------------------------------------*/
static void ShowVector(VECTOR *pVSrc)
{
WORD i;
    for(i = 0; i < pVSrc->size; i++ )
        printf(" %04x", pVSrc->data[i]);
    printf("\n");
}

#if EEUCLID
/*----------------------------------------------------------------------*/
/*      ShowEuclid                                                      */
/*----------------------------------------------------------------------*/
static void ShowEuclid(EUCLID *pESrc)
{
WORD i;
    for(i = 0; i < pESrc->indx; i++){
        printf(" %04x", pESrc->data[i]);}
    printf("|");
    for( ; i < pESrc->size; i++){
        printf("%04x ", pESrc->data[i]);}
    printf("\n");
}
#endif

