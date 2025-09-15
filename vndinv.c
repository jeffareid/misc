/*----------------------------------------------------------------------*/
/*      vndinv.c    test vandermonde inversion                          */
/*                                                                      */
/*      Jeff Reid   2025SEP14  17:30                                    */
/*----------------------------------------------------------------------*/
/*      equates                                                         */
/*----------------------------------------------------------------------*/
#define _CRT_SECURE_NO_WARNINGS 1       /* disable sscanf warnings */

#include <stdio.h>
#include <string.h>
#include <memory.h>

typedef unsigned  char  BYTE;

/*                                      ** display matrix inv */
#define DISPLAYI 0

typedef struct{                         /* matrix structure */
    int  nrows;
    int  ncols;
    BYTE  data[256*256*2];
}MATRIX;

/*----------------------------------------------------------------------*/
/*      data                                                            */
/*----------------------------------------------------------------------*/
static int aiGA[60] =                   /* GF's and Alpha's */
   {0x11b,0x03,0x11d,0x02,0x12b,0x02,0x12d,0x02,
    0x139,0x03,0x13f,0x03,0x14d,0x02,0x15f,0x02,
    0x163,0x02,0x165,0x02,0x169,0x02,0x171,0x02,
    0x177,0x03,0x17b,0x0b,0x187,0x02,0x18b,0x0b,
    0x18d,0x02,0x19f,0x03,0x1a3,0x03,0x1a9,0x02,
    0x1b1,0x07,0x1bd,0x07,0x1c3,0x02,0x1cf,0x02,
    0x1d7,0x07,0x1dd,0x07,0x1e7,0x02,0x1f3,0x0d,
    0x1f5,0x02,0x1f9,0x03};

static BYTE     abExp[512];             /* antilog table */
static BYTE     abLog[256];             /* log table */

static BYTE     abId[516];              /* array for MatrixInv */

static int      iGF;                    /* Galios Field Polynomial */
static BYTE     bAlpha;                 /* Alpha for this field */

/*----------------------------------------------------------------------*/
/*      code                                                            */
/*----------------------------------------------------------------------*/
static int      MatrixInv(MATRIX *, MATRIX *);
static void     MatrixMpy(MATRIX *, MATRIX *, MATRIX *);
static BYTE     GFAdd(BYTE, BYTE);
static BYTE     GFSub(BYTE, BYTE);
static BYTE     GFMpy(BYTE, BYTE);
static BYTE     GFDiv(BYTE, BYTE);
static BYTE     GFPow(BYTE, BYTE);
static void     InitGF(void);
static void     DoUser(void);


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
    abId[256] =  1;

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
        memcpy(p1,  &abId[256-i], iNCol*sizeof(BYTE)); /* add ID  part */
        p1  += iNCol;}

/*      normalize according to left side */
/*      results in inverse matrix in right size */

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
            ;}


/*      divide row to produce a one */

        p0  = pMDst->data+j*iNCol2;     /* p0 = ptr to  start of row */
        for(p2  = p0+iNCol2; p0 != p2; p0++){
            *p0 = GFDiv(*p0, bMod);}

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
    return(0);
}

/*----------------------------------------------------------------------*/
/*      MatrixMpy(pMDst, pMSrc0, pmSrc1) matrix multiply                */
/*----------------------------------------------------------------------*/
static void MatrixMpy(MATRIX *pMDst, MATRIX *pMSrc0, MATRIX *pMSrc1)
{
int i, j, k;
int drows, dcols, inner;
BYTE *pucd;                             /* ptr to dst */
BYTE *pucs0, *pucs1;                    /* ptr to src */

    drows = pMSrc0->nrows;              /* calc dest params */
    dcols = pMSrc1->ncols;
    inner = pMSrc0->ncols;              /* inner product size */

    pMDst->nrows = drows;               /* init dest */
    pMDst->ncols = dcols;
    memset(pMDst->data, 0, drows*dcols);

    pucd = pMDst->data;                 /* do the mpy */
    for(k = 0; k < drows; k++){
        for(j = 0; j < dcols; j++){
            pucs0 = &pMSrc0->data[k*inner];
            pucs1 = &pMSrc1->data[j];
            for(i = 0; i < inner; i++){
                *pucd ^= GFMpy(*pucs0, *pucs1);
                pucs0 += 1;
                pucs1 += dcols;}
            pucd += 1;}}
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
        return(0);}
    if(b0 == 0)
        return(0);
    return(abExp[255+(int)abLog[b0]-(int)abLog[b1]]);
}

/*----------------------------------------------------------------------*/
/*      GFPow(b0, b1)           b0^b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFPow(BYTE b0, BYTE b1)
{
BYTE b;
    b = 1;
    while(b1){
        if(b1&1)
            b = GFMpy(b, b0);
        b0 = GFMpy(b0, b0);
        b1 >>= 1;}
    return(b);
}

/*----------------------------------------------------------------------*/
/*      GFMpy0(b0,b1)           b0*b1       using low level math        */
/*----------------------------------------------------------------------*/
static BYTE GFMpy0(BYTE b0, BYTE b1)
{
int i;
int product;
    product = 0;
    for(i = 0; i < 8; i++){
        product <<= 1;
        if(product & 0x100){
            product ^= iGF;}
        if(b0 & 0x80u){
            product ^= b1;}
        b0 <<= 1;}
    return((BYTE)product);
}

/*----------------------------------------------------------------------*/
/*      InitGF  Initialize Galios Stuff                                 */
/*----------------------------------------------------------------------*/
static void InitGF(void)
{
BYTE b;
int i;

    b = 1;
    for(i = 0; i < 512; i++){           /* init abExp[] */
        abExp[i] = b;
        b = GFMpy0(b, bAlpha);}

    abLog[0] = 0xff;                    /* init abLog[] */
    for(i = 0; i < 255; i++){
        abLog[abExp[i]] = i;}
}

/*----------------------------------------------------------------------*/
/*      main                                                            */
/*----------------------------------------------------------------------*/
int main()
{
int i, j, k;
MATRIX m1, m2, m3;

/* select Galios Field */
    iGF    = aiGA[2];
    bAlpha = (BYTE)(aiGA[3]);
    InitGF();
    for(k = 1; k < 256; k++){
        m1.nrows = k;
        m1.ncols = k;
        for(j = 0; j < m1.nrows; j++)
            for(i = 0; i < m1.ncols; i++)
                m1.data[(j*m1.ncols)+i] = GFPow(j, i);
        MatrixInv(&m2, &m1);        /* invert Vandermonde matrix */
        MatrixMpy(&m3, &m2, &m1);   /* check it was inverted */
        for(i = 0; i < k; i++) {
            if(m3.data[k*i+i] != 1){
                goto err0;
            }
            m3.data[k*i+i] = 0;
        }
        for(j = 0; j < k; j++){
            for(i = 0; i < k; i++){
                if(m3.data[k*j+i] != 0){
                    goto err0;
                }
            }
        }
        printf(".");
    }
    printf("\npass\n");
    return 0;
err0:
    printf("\nfail\n");
    return 0;
}
