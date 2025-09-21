/*----------------------------------------------------------------------*/
/*      bch63.c         63 bit codeword, GF(2^6)                        */
/*                                                                      */
/*      Jeff Reid       2021JAN09 11:30                                 */
/*----------------------------------------------------------------------*/
#include <intrin.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned char      BYTE;
typedef unsigned short     WORD;
typedef unsigned int       DWORD;
typedef unsigned long long QWORD;

/*                              ** GF(2^6): x^6 + x + 1 */
#define POLY (0x43)

/*                              ** GF(2^6) primitive */
#define ALPHA (0x2)

/*                              ** 4 error correction, 8 syndromes */
#define NSYN (8)

/*                              ** use extended euclid algorithm */
#define EEUCLID 0

/*                              ** display euclid stuff */
#define DISPLAYE 0

/*                              ** display error locator poly */
#define DISPLAYP 1

typedef struct{                 /* vector structure */
    BYTE  size;
    BYTE  data[15];
}VECTOR;

#if EEUCLID
typedef struct{                 /* euclid structure */
    BYTE  size;                 /* # of data bytes */
    BYTE  indx;                 /* index to right side */
    BYTE  data[NSYN+2];         /* left and right side data */
}EUCLID;
#endif

static __m128i poly;            /* generator poly */
static __m128i invpoly;         /* 2^64 / POLY */

static BYTE log2[64];           /* log2 table */
static BYTE alog2[64];          /* alog2 table */
static BYTE minplyf[128];       /* minimum polynomial flags */
static BYTE minply[64];         /* minimum polynomials */
static BYTE mincnt;             /* # of minimum polymials */

static __m128i psyn[NSYN];      /* syndrome polynomials */
static __m128i invsyn[NSYN];    /* inverse syndrome polynomials */

static __m128i msg;             /* msg.m128i_u64[0] is test message */
static __m128i par;             /* par.m128i_u64[1] is parities */

static QWORD   encmsg;          /* copy of encoded message */

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

static void Tbli(void);
static void GenPoly(void);
static void Encode(void);
static void GenSyndromes(void);
static void GenpErrors(void);
static void GenOffsets(void);
static void FixErrors(void);
static int  Poly2Root(VECTOR *, VECTOR *);
static BYTE GFPwr(BYTE, BYTE);
static BYTE GFMpy(BYTE, BYTE);
static BYTE GFDiv(BYTE, BYTE);
static void ShowVector(VECTOR *);
#if EEUCLID
static void ShowEuclid(EUCLID *);
#endif
static void InitCombination(int[], int, int);
static int  NextCombination(int[], int, int);

#define GFAdd(a, b) (a^b)
#define GFSub(a, b) (a^b)

/*----------------------------------------------------------------------*/
/*      main                                                            */
/*----------------------------------------------------------------------*/
main()
{
int ptn[4];                     /* error bit indexes */
int n;                          /* number of errors to test */
int i;
    Tbli();                     /* init tables */
    GenPoly();                  /* generate poly info */
    msg.m128i_u64[0] = 0x123456789a000000ull;       /* test message */
    Encode();                   /* encode message */
    encmsg = msg.m128i_u64[0];
    for(n = 1; n <= 4; n++){    /* test 1 to 4 bit error patterns */
        InitCombination(ptn, n, 63);
        while(NextCombination(ptn, n, 63)){
            for(i = 0; i < n; i++)
                msg.m128i_u64[0] ^= 1ull<<ptn[i];
            GenSyndromes();     /* generate syndromes */
            GenpErrors();       /* generate error location info */
            GenOffsets();       /* convert to offsets */
            FixErrors();        /* correct error bits */
            if(encmsg != msg.m128i_u64[0]){
                printf("failed\n");
                return 0;
            }
        }
    }
    printf("passed\n");
    return 0;
}

/*----------------------------------------------------------------------*/
/*      Tbli                                                            */
/*----------------------------------------------------------------------*/
static void Tbli()              /* init tables */
{
BYTE *p0, *p1;
int d0;
    d0 = 01;                        /* init alog2 table */
    p0 = alog2;
    for(p1 = p0+64; p0 < p1;){
        *p0++ = d0;
        d0 <<= 1;
        if(d0 & 0x40)
            d0 ^= POLY;}

    p0 = alog2;                     /* init log2 table */
    p1 = log2;
    *p1 = 0;
    for(d0 = 0; d0 < 63; d0 += 1)
        *(p1+*p0++) = d0;
}

/*----------------------------------------------------------------------*/
/*      GenPoly                                                         */
/*----------------------------------------------------------------------*/
static void GenPoly(void)
{
QWORD M;
QWORD N;
QWORD Q;
DWORD x;
BYTE mpoly;                     /* test polynomial */
BYTE sum;                       /* sum, looking for zeroes */
BYTE apwr;                      /* alpha to power */
BYTE i,j;

    /* find minimum and non-duplicate polynomials for m[1] -> m[NSYN-1] */
    for(i = 1; i < NSYN; i++){
        apwr = GFPwr(ALPHA,i);
        for(mpoly = 0x02; mpoly <= 0x7f ; mpoly++){
            sum = 0;
            for(j = 0; j <= 6; j++){
                if(mpoly&(1<<j))
                    sum ^= GFPwr(apwr,j);
            }
            if(sum == 0){
                if(minplyf[mpoly] != 0)
                    continue;
                minplyf[mpoly] = 1;
                minply[mincnt] = mpoly;
                mincnt += 1;
                break;
            }
        }
    }

    poly.m128i_u64[0] = minply[0];
    for(i = 1; i < mincnt; i++){
        poly.m128i_u64[1] = minply[i];
        poly = _mm_clmulepi64_si128(poly, poly, 0x01);
    }

    /* generate inverse of encoding polynomial */
    _BitScanReverse64(&x, poly.m128i_u64[0]);
    M = 1ull<<x;
    N = M>>1;
    Q = 0x0ull;
    for(i = 0; i < 65-x; i++){
        N <<= 1;
        Q <<= 1;
        if(N&M){
            Q |= 1;
            N ^= (poly.m128i_u64[0]);
        }
    }
    invpoly.m128i_u64[0] = Q;
}

/*----------------------------------------------------------------------*/
/*      Encode                                                          */
/*----------------------------------------------------------------------*/
static void Encode(void)
{
    par = _mm_clmulepi64_si128(msg, invpoly, 0x00); /* par[1] = quotient */
    par = _mm_clmulepi64_si128(par, poly, 0x01);    /* par[0] = product */
    par.m128i_u64[0] ^= msg.m128i_u64[0];           /* par[0] = remainder */
    msg.m128i_u64[0] |= par.m128i_u64[0];           /* msg[0] = encoded message */
}

/*----------------------------------------------------------------------*/
/*      GenSyndromes                                                    */
/*----------------------------------------------------------------------*/
static void GenSyndromes(void)
{
QWORD M;
DWORD x;
BYTE i, ap, s;
    vSyndromes.size = NSYN;
    for(i = 0; i < NSYN; i++){
        M = msg.m128i_u64[0];
        ap = GFPwr(ALPHA, i+1);
        s = 0;
        while(M){
            _BitScanReverse64(&x, M);
            M ^= 1ull<<x;
            s ^= GFPwr(ap, x);
        }
        vSyndromes.data[i] = s;
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
BYTE    bME;                            /* max errors possible */
BYTE    bQuot;                          /* quotient */

/*      E0 = initial ED: E0.R[-1] = x^MAXERR, E0.A[0] = 1 */
    E0.size = vSyndromes.size+2;
    E0.indx = vSyndromes.size+1;
    E0.data[0] = 1;
    memset(&E0.data[1], 0, vSyndromes.size*sizeof(BYTE));
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
            memcpy(&pER->data[0], &pER->data[1], (pER->size-1)*sizeof(BYTE));
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
            memcpy(&pED->data[0], &pED->data[1], (pED->size-1)*sizeof(BYTE));
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
BYTE i, j, n;
BYTE L, m;
BYTE b, d;
BYTE db;

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
            continue;
        }
        d = vSyndromes.data[n];         /* calculate discrepancy */
        for(i = 1; i <= L; i++){
            d = GFAdd(d, GFMpy(vC.data[(vC.size - 1)- i], vSyndromes.data[n-i]));}
        if(d == 0){                     /* if 0 increment m, continue */
            m += 1;
            continue;}
        vT.size = vC.size;              /* vT = vC */
        memcpy(vT.data, vC.data, vC.size);
        db = GFDiv(d,b);                /* db = (d/b) */
        vBx.size = vB.size+m;           /* Bx = x^m B */
        memcpy(vBx.data, vB.data, vB.size);
        memset(&vBx.data[vB.size], 0, m);
        for(i = 0; i < vBx.size; i++){  /* Bx *= db */
            vBx.data[i] = GFMpy(vBx.data[i], db);}
        j = vBx.size - vC.size;         /* right shift vBx or vC */
        if(((char)j) > 0){
            for(i = vBx.size; i > j; ){
                i--;
                vC.data[i] = vC.data[i-j];}
            memset(vC.data, 0, j);
            vC.size += j;}
        else if(((char)j) < 0){
            j = -j;
            for(i = vC.size; i > j; ){
                i--;
                vBx.data[i] = vBx.data[i-j];}
            memset(vBx.data, 0, j);
            vBx.size += j;}
        for(i = 0; i < vC.size; i++){   /* C -= Bx */
            vC.data[i] = GFSub(vC.data[i], vBx.data[i]);}
        if(n < 2*L){                    /* if L not increasing */
            m += 1;
            continue;}
        vB.size = vT.size;              /*   B = T */
        memcpy(vB.data, vT.data, vT.size);
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
BYTE i;
    vOffsets.size = vLocators.size;
    for(i = 0; i < vLocators.size; i++){
        vOffsets.data[i] = log2[vLocators.data[i]];
    }
}

/*----------------------------------------------------------------------*/
/*      FixErrors                                                       */
/*----------------------------------------------------------------------*/
static void FixErrors()
{
BYTE i;
    for(i = 0; i < vOffsets.size; i++)
        msg.m128i_u64[0] ^= 1ull<<vOffsets.data[i];
}

/*----------------------------------------------------------------------*/
/*      Poly2Root(pVDst, pPSrc)         find roots of poly              */
/*----------------------------------------------------------------------*/
static int Poly2Root(VECTOR *pVDst, VECTOR *pPSrc)
{
BYTE    bLcr;                           /* current locator */
BYTE    bSum;                           /* current sum */
BYTE    bV;                             /* index to pVDst */
BYTE    i,j;

    pVDst->size = pPSrc->size-1;        /* set dest size */

    if(!pVDst->size)                    /* exit if null */
        return(0);

    bV   = 0;
    bLcr = 1;
    for(j = 0; j < 63;  j++){
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
/*      GFPwr                                                           */
/*----------------------------------------------------------------------*/
static BYTE GFPwr(BYTE m0, BYTE m1)
{
    return alog2[(BYTE)((log2[m0]*(DWORD)m1)%63)];
}

/*----------------------------------------------------------------------*/
/*      GFMpy                                                           */
/*----------------------------------------------------------------------*/
static BYTE GFMpy(BYTE m0, BYTE m1) /* multiply */
{
int m2;
    if(0 == m0 || 0 == m1)
        return(0);
    m2 = log2[m0] + log2[m1];
    if(m2 > 63)
        m2 -= 63;
    return(alog2[m2]);
}

/*----------------------------------------------------------------------*/
/*      GFDiv                                                           */
/*----------------------------------------------------------------------*/
static BYTE GFDiv(BYTE m0, BYTE m1) /* divide */
{
int m2;
    if(0 == m0)
        return(0);
    m2 = log2[m0] - log2[m1];
    if(m2 < 0)
        m2 += 63;
    return(alog2[m2]);
}

/*----------------------------------------------------------------------*/
/*      ShowVector                                                      */
/*----------------------------------------------------------------------*/
static void ShowVector(VECTOR *pVSrc)
{
BYTE    i;
    for(i = 0; i < pVSrc->size; ){
        printf(" %02x", pVSrc->data[i]);
        i++;
        if(0 == (i&0xf)){
            printf("\n");}}
    printf("\n");
}


#if EEUCLID
/*----------------------------------------------------------------------*/
/*      ShowEuclid                                                      */
/*----------------------------------------------------------------------*/
static void ShowEuclid(EUCLID *pESrc)
{
BYTE    i;
    for(i = 0; i < pESrc->indx; i++){
        printf(" %02x", pESrc->data[i]);}
    printf("|");
    for( ; i < pESrc->size; i++){
        printf("%02x ", pESrc->data[i]);}
    printf("\n");
}
#endif

/*----------------------------------------------------------------------*/
/*      InitCombination - init combination                              */
/*----------------------------------------------------------------------*/
void InitCombination(int a[], int k, int n) {
    for(int i = 0; i < k; i++)
        a[i] = i;
    --a[k-1];
}

/*----------------------------------------------------------------------*/
/*      NextCombination - generate next combination                     */
/*----------------------------------------------------------------------*/
int NextCombination(int a[], int k, int n) {
int pivot = k - 1;
    while (pivot >= 0 && a[pivot] == n - k + pivot)
        --pivot;
    if (pivot == -1)
        return 0;
    ++a[pivot];
    for (int i = pivot + 1; i < k; ++i)
        a[i] = a[pivot] + i - pivot;
    return 1;
}