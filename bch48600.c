/*----------------------------------------------------------------------*/
/*      bch48600.c      bch((48600,48408)  DVBS2 code                   */
/*                                                                      */
/*      Jeff Reid       2021JAN09 11:00                                 */
/*----------------------------------------------------------------------*/
#include <intrin.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned char      BYTE;
typedef unsigned short     WORD;
typedef unsigned long      DWORD;
typedef unsigned long long QWORD;

/* GF(2^16) = x^16 + x^5 + x^3 + x^2 + 1 */
#define POLY 0x1002D

/* GF(2^16) primitive */
#define ALPHA 0x0002

/* BCH(48600, 48408) */
#define BN 48600
#define BK 48408
#define BP 192

/* number of suyndromes */
#define NSYN 24

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
static WORD gflog2[65536];      /* gflog2 table */
static WORD gfexp2[65536];      /* gfexp2 table */

static DWORD minply[12];        /* minimum polynomials */
static WORD mincnt;             /* # of minimum polymials */
static BYTE minplyf[65536];     /* minimum polynomial flags */

static QWORD polytbl[256][3];   /* encode poly table */

static WORD syntbl[256][NSYN];  /* syndrome table */

static BYTE msg[(BN+7)/8];      /* encoded message */

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
    GFInit();                   /* init GF(16) tables */
    GenMinPoly();               /* generate minimum polynomials */
    GenPolyTbl();               /* generate 192 bit poly table */
    GenSynTbl();                /* generate syndrome table */
    for(i = 0; i < (BK/8); i++) /* gnerate test message */
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
    msg[   0] ^= 0x80;          /* test 12 error bit case */
    msg[ 506] ^= 0x40;
    msg[1012] ^= 0x20;
    msg[1518] ^= 0x10;
    msg[2024] ^= 0x08;
    msg[2530] ^= 0x04;
    msg[3543] ^= 0x02;
    msg[4050] ^= 0x01;
    msg[4556] ^= 0x08;
    msg[5062] ^= 0x04;
    msg[5568] ^= 0x02;
    msg[6074] ^= 0x01;
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
/*      GenMinPoly      generate 12 mininum polynomials                 */
/*----------------------------------------------------------------------*/
static void GenMinPoly(void)
{
DWORD poly;                     /* test polynomial */
DWORD sum;                      /* sum, looking for zeroes */
WORD apwr;                      /* alpha to power */
WORD i,j;

    /* find 12 minimum polynomials for 24 powers of 2 */
    /* result matches the gif polynomials */
    i = 0;
    do{
        apwr = GFPwr(2,++i);
        for(poly = 0x10001; poly <= 0x1ffff ; poly++){
            sum = 0;
            for(j = 0; j <= 16; j++){
                if(poly&(1<<j))
                    sum ^= GFPwr(apwr,j);
            }
            if(sum == 0){
                if(minplyf[poly-0x10000] != 0)
                    continue;
                minplyf[poly-0x10000] = 1;
                minply[mincnt++] = poly;
                break;
            }
        }
    }while(i < NSYN);
}

/*----------------------------------------------------------------------*/
/*      GenPolyTbl     generate 192 bit polynomial table                */
/*----------------------------------------------------------------------*/
static void GenPolyTbl(void)
{
__m128i p[4];                   /* sub-products and product */
QWORD poly[3];                  /* 192 bit polynomial */
QWORD t[3];                     /* for table generation */
QWORD b, q, i;                  /* byte, quotient bit, i */

/*  gnerate two 49 bit products */
    p[0].m128i_u64[0] = minply[ 0];
    p[0].m128i_u64[1] = minply[ 1];
    p[0]              = _mm_clmulepi64_si128(p[0], p[0], 0x01);
    p[0].m128i_u64[1] = minply[ 2];
    p[0]              = _mm_clmulepi64_si128(p[0], p[0], 0x01);
    p[1].m128i_u64[0] = minply[ 3];
    p[1].m128i_u64[1] = minply[ 4];
    p[1]              = _mm_clmulepi64_si128(p[1], p[1], 0x01);
    p[1].m128i_u64[1] = minply[ 5];
    p[1]              = _mm_clmulepi64_si128(p[1], p[1], 0x01);
/* generate 97 bit product */   
    p[2]              =  _mm_clmulepi64_si128(p[0], p[1], 0x00);
/*  gnerate two 49 bit products */
    p[0].m128i_u64[0] = minply[ 6];
    p[0].m128i_u64[1] = minply[ 7];
    p[0]              = _mm_clmulepi64_si128(p[0], p[0], 0x01);
    p[0].m128i_u64[1] = minply[ 8];
    p[0]              = _mm_clmulepi64_si128(p[0], p[0], 0x01);
    p[1].m128i_u64[0] = minply[ 9];
    p[1].m128i_u64[1] = minply[10];
    p[1]              = _mm_clmulepi64_si128(p[1], p[1], 0x01);
    p[1].m128i_u64[1] = minply[11];
    p[1]              = _mm_clmulepi64_si128(p[1], p[1], 0x01);
/* generate 97 bit product */   
    p[3]              = _mm_clmulepi64_si128(p[0], p[1], 0x00);
/* generate lower 128 bit product */
    p[0]              = _mm_clmulepi64_si128(p[2], p[3], 0x00);
    poly[2]       = p[0].m128i_u64[0];
    poly[1]       = p[0].m128i_u64[1];
/* generate middle 128 bit products */
    p[0]              = _mm_clmulepi64_si128(p[2], p[3], 0x01);
    poly[1]      ^= p[0].m128i_u64[0];
    poly[0]       = p[0].m128i_u64[1];
    p[0]              = _mm_clmulepi64_si128(p[2], p[3], 0x10);
    poly[1]      ^= p[0].m128i_u64[0];
    poly[0]      ^= p[0].m128i_u64[1];
/* generate upper 128 bit product */
    p[0]              = _mm_clmulepi64_si128(p[2], p[3], 0x11);
    poly[0]      ^= p[0].m128i_u64[0];
/* generate polytbl */  
    for(b = 0x000; b < 0x100; b++){
        t[0] = b<<56;
        t[1] = 0;
        t[2] = 0;
        for(i = 0; i < 8; i++){
            q = t[0]>>63;
            t[0] = (t[0]<<1)|(t[1]>>63);
            t[1] = (t[1]<<1)|(t[2]>>63);
            t[2] = (t[2]<<1);
            if(q != 0){
                t[0] ^= poly[0];
                t[1] ^= poly[1];
                t[2] ^= poly[2];
            }
        }
        polytbl[b][0] = t[0];
        polytbl[b][1] = t[1];
        polytbl[b][2] = t[2];
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
QWORD r[3];                     /* remainder */
QWORD i, j;
    r[2] = r[1] = r[0] = 0;
    for(j = 0; j < BK/8; j += 1){
        i = (r[0]>>56)^msg[j];  /* table index */
        r[0] = (r[0]<<8)|(r[1]>>56);
        r[1] = (r[1]<<8)|(r[2]>>56);
        r[2] = (r[2]<<8);
        r[0] ^= polytbl[i][0];
        r[1] ^= polytbl[i][1];
        r[2] ^= polytbl[i][2];
    }
    msg[(BK/8)+ 0] = (BYTE)(r[0]>>56);
    msg[(BK/8)+ 1] = (BYTE)(r[0]>>48);
    msg[(BK/8)+ 2] = (BYTE)(r[0]>>40);
    msg[(BK/8)+ 3] = (BYTE)(r[0]>>32);
    msg[(BK/8)+ 4] = (BYTE)(r[0]>>24);
    msg[(BK/8)+ 5] = (BYTE)(r[0]>>16);
    msg[(BK/8)+ 6] = (BYTE)(r[0]>> 8);
    msg[(BK/8)+ 7] = (BYTE)(r[0]>> 0);
    msg[(BK/8)+ 8] = (BYTE)(r[1]>>56);
    msg[(BK/8)+ 9] = (BYTE)(r[1]>>48);
    msg[(BK/8)+10] = (BYTE)(r[1]>>40);
    msg[(BK/8)+11] = (BYTE)(r[1]>>32);
    msg[(BK/8)+12] = (BYTE)(r[1]>>24);
    msg[(BK/8)+13] = (BYTE)(r[1]>>16);
    msg[(BK/8)+14] = (BYTE)(r[1]>> 8);
    msg[(BK/8)+15] = (BYTE)(r[1]>> 0);
    msg[(BK/8)+16] = (BYTE)(r[2]>>56);
    msg[(BK/8)+17] = (BYTE)(r[2]>>48);
    msg[(BK/8)+18] = (BYTE)(r[2]>>40);
    msg[(BK/8)+19] = (BYTE)(r[2]>>32);
    msg[(BK/8)+20] = (BYTE)(r[2]>>24);
    msg[(BK/8)+21] = (BYTE)(r[2]>>16);
    msg[(BK/8)+22] = (BYTE)(r[2]>> 8);
    msg[(BK/8)+23] = (BYTE)(r[2]>> 0);
}

/*----------------------------------------------------------------------*/
/*      GenSyndromes                                                    */
/*----------------------------------------------------------------------*/
static void GenSyndromes(void)
{
WORD i, j, apwr;
    vSyndromes.size = NSYN;
    memset(vSyndromes.data, 0, NSYN*sizeof(WORD));
    for(j = 0; j < BN/8; j++){
        for(i = 0; i < NSYN; i++){
            apwr = GFPwr(2, (i+1)<<3);
            vSyndromes.data[i] = GFMpy(vSyndromes.data[i], apwr) ^
                                 syntbl[msg[j]][i];
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
        if(((char)j) > 0){
            for(i = vBx.size; i > j; ){
                i--;
                vC.data[i] = vC.data[i-j];}
            memset(vC.data, 0, j*sizeof(WORD));
            vC.size += j;}
        else if(((char)j) < 0){
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
        vOffsets.data[i] = BN-1-gflog2[vLocators.data[i]];
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
/*      GFInit          init GF(16) tables                              */
/*----------------------------------------------------------------------*/
static void GFInit(void)        /* init tables */
{
WORD *p0, *p1;
int d0;
d0 = 01;                        /* init gfexp2 table */
p0 = gfexp2;
for(p1 = p0+65536; p0 < p1;){
  *p0++ = d0;
  d0 <<= 1;
  if(d0 & 0x10000)
    d0 ^= POLY;}

p0 = gfexp2;                        /* init gflog2 table */
p1 = gflog2;
*p1 = 0;
for(d0 = 0; d0 < 65535; d0 += 1)
  *(p1+*p0++) = d0;
}

/*----------------------------------------------------------------------*/
/*      GFPwr - exponentiate                                            */
/*----------------------------------------------------------------------*/
static WORD GFPwr(WORD m0, WORD m1)
{
    return gfexp2[(WORD)((gflog2[m0]*(DWORD)m1)%65535)];
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
    if(m2 > 65535)
        m2 -= 65535;
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
        m2 += 65535;
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

