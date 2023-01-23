/*----------------------------------------------------------------------*/
/*      m8to44f.c       find mapping (16^2) field given (2^8) field     */
/*                                                                      */
/*      Jeff Reid       11MAR27  20:00                                  */
/*----------------------------------------------------------------------*/
#include <memory.h>
#include <stdio.h>

typedef unsigned char BYTE;

/*                                      ** only allow alpha = 1x + 0 */
#define CMAX 0x01
#define DMAX 0x00

typedef struct{                         /* polynomial structure */
    BYTE  size;
    BYTE  data[11];
}POLY;

/*----------------------------------------------------------------------*/
/*      data                                                            */
/*----------------------------------------------------------------------*/
static BYTE log4[16];                   /* 4-bit log table */
static BYTE exp4[16];                   /* 4-bit exp table */
static BYTE prm44[256];                 /* prime gf(16) poly table */

static int  tgf8[30] =                  /* GF(2^8) polynomials */
   {0x11b,0x11d,0x12b,0x12d,0x139,0x13f,0x14d,0x15f,
    0x163,0x165,0x169,0x171,0x177,0x17b,0x187,0x18b,
    0x18d,0x19f,0x1a3,0x1a9,0x1b1,0x1bd,0x1c3,0x1cf,
    0x1d7,0x1dd,0x1e7,0x1f3,0x1f5,0x1f9};

static BYTE tal8[30] =                  /* GF(2^8) alphas */
   {0x003,0x002,0x002,0x002,0x003,0x003,0x002,0x002,
    0x002,0x002,0x002,0x002,0x003,0x00b,0x002,0x00b,
    0x002,0x003,0x003,0x002,0x007,0x007,0x002,0x002,
    0x007,0x007,0x002,0x00d,0x002,0x003};

static int  t8idx;                      /* table index */
static int  gf8;                        /* GF( 2^8) polynomial */
static BYTE al8;                        /* GF( 2^8) alpha */
static BYTE bpad0;                      /* pad (unused) */

static BYTE tgf4[3]={0x13, 0x19, 0x1f}; /* GF(2^4) polynomials */
static BYTE tal4[3]={0x02, 0x02, 0x03}; /* GF(2^4) alphas */
static BYTE gf4;
static BYTE al4;
static int  t4idx;                      /* table index */

static BYTE a;                          /* hi nibble of test value */
static BYTE b;                          /* lo nibble of test value */
static BYTE c;                          /* hi nibble of test alpha */
static BYTE d;                          /* lo nibble of test alpha */
static BYTE q;                          /* quotient nibble */
static BYTE gf0, gf1;                   /* poly nibbles */

static BYTE exp08[256];                 /* GF( 2^8) exp[] table */
static BYTE exp44[256];                 /* GF(16^2) exp[] table to test */

static POLY p44[8];                     /* compatible GF(16^2) polynomails alpha == 1x + 0 */
static POLY p4x[8];                     /* compatible GF(16^2) polynomails alpha != 1x + 0 */
static POLY p1, p2;                     /* used for polympy() */
static int  p44idx;                     /* index into p44 */
static int  p4xidx;                     /* index into p4x */

/*----------------------------------------------------------------------*/
/*      code                                                            */
/*----------------------------------------------------------------------*/

static int  TestMap(void);
static void InitTables(void);
static void PolyMpy(POLY *, POLY *, POLY *);
static BYTE Mpy4(BYTE, BYTE);
static BYTE Div4(BYTE, BYTE);
static BYTE Mp4x(BYTE, BYTE);
static BYTE Mp8x(BYTE, BYTE);

/*----------------------------------------------------------------------*/
/*      main                                                            */
/*----------------------------------------------------------------------*/
main(int argc, char **argv)
{
int i, j, k;

    for(t4idx = 0; t4idx < sizeof(tal4); t4idx++){  /* select 4 bit field */
        gf4 = tgf4[t4idx];
        al4 = tal4[t4idx];

        for(t8idx = 0; t8idx < sizeof(tal8); t8idx++){ /* select 8 bit field */
            gf8 = tgf8[t8idx];
            al8 = tal8[t8idx];
            if(al8 != 0x02){                        /* if al8 != 0x02, not mappable */
                continue;}

            p44idx = 0;                             /* reset index into p44 */
            p4xidx = 0;

            InitTables();                           /* init tables */

            printf("gf4 %02x  gf8 %03x", gf4, gf8);

            for(i = 0; i < 0x100; i++){             /* init prime table */
                prm44[i] = i;}

            for(b = 0; b < 0x10; b++){              /* zero out non-primes polynomials */
                for(a = 0; a < 0x10; a++){
                    i = ((a^b)<<4)^(Mpy4(a,b));     /*  i = (x - a)(x - b) */
                    prm44[i] = 0;}}

            for(k = 0; k < 0x100; k++){             /* for all polynomials */
                if(!prm44[k])continue;              /*  skip if non-prime polynomial */
                gf1 = (prm44[k]>>4);                /* poly = 1 x^2 + gf1 x + gf0 */
                gf0 = (prm44[k]&0xf);
                for(c = 1; c <= CMAX; c++){         /* for each alpha (cx + d) */
                    for(d = 0; d <= DMAX; d++){
                      i = 0x01;                     /*  i = alpha^0 */
                      for(j = 0; j < 0x100; j++){   /*  cycle through powers of alpha */
                          a = (i>>4);               /*  i *= alpha */
                          b = (i&0xf);
                          i = ((Mpy4(Mpy4(a,c),gf1)^Mpy4(a,d)^Mpy4(b,c))<<4)^Mpy4(Mpy4(a,c),gf0)^Mpy4(b,d);
                          if(i == 0x01)break;}      /*  break if cycled back to alpha^0 */
                      if(j == 0xfe){                /* if this alpha is the primitive, */
                          TestMap();                /*   test for compatability */
                          goto nextpoly;}}}
nextpoly:       continue;}
            if(p4xidx || p44idx){
                printf(" gf44 ");}
            if(p4xidx){
                for(i = 0; i < p4xidx; i++){
                    printf("1%1x%1x:%1x%1x ",
                        p4x[i].data[1], p4x[i].data[2],
                        p4x[i].data[3], p4x[i].data[4]);}
                printf(" ");}
            if(p44idx){
                p1.size = 1;
                p1.data[0] = 1;
                for(i = 0; i < p44idx; i++){
                    printf("1%1x%1x ",
                        p44[i].data[1], p44[i].data[2]);
                    PolyMpy(&p2, &p1, &p44[i]);
                    p1 = p2;}
                printf(" prod ");
                j = 0;
                for(i = 0; i < p1.size; i++){
                    printf("%1x", p1.data[i]);
                    j <<= 1;
                    j ^= p1.data[i];}
                printf(" %03x", j);}
            printf("\n");}
        printf("\n");}
    return(0);
}

/*----------------------------------------------------------------------*/
/*      TestMap                                                         */
/*                                                                      */
/*      if (1 x^2 + gf1 x + gf0) compatable with gf8, then              */
/*      The first 8 entries of exp44[] map GF(2^8)->GF(16^2)            */
/*      define [i.b] to mean exp44[i] bit b:                            */
/*      define [c] to mean GF(2^8) byte bit c:                          */
/*      then conversion GF(2^8)->GF(16^2) is this matrix multiply:      */
/*                                                                      */
/*      [7.7] [6.7] [5.7] [4.7] [3.7] [2.7] [1.7] [0.7]   [7]           */
/*      [7.6] [6.6] [5.6] [4.6] [3.6] [2.6] [1.6] [0.6]   [6]           */
/*      [7.5] [6.5] [5.5] [4.5] [3.5] [2.5] [1.5] [0.5]   [5]           */
/*      [7.4] [6.4] [5.4] [4.4] [3.4] [2.4] [1.4] [0.4]   [4]           */
/*      [7.3] [6.3] [5.3] [4.3] [3.3] [2.3] [1.3] [0.3]   [3]           */
/*      [7.2] [6.2] [5.2] [4.2] [3.2] [2.2] [1.2] [0.2]   [2]           */
/*      [7.1] [6.1] [5.1] [4.1] [3.1] [2.1] [1.1] [0.1]   [1]           */
/*      [7.0] [6.0] [5.0] [4.0] [3.0] [2.0] [1.0] [0.0]   [0]           */
/*----------------------------------------------------------------------*/
static int TestMap(void)
{
int i, j, k;
BYTE *pb;
BYTE b44;

    pb = exp44;                     /* generate exp table for this poly */
    i = 0x01;
    *pb++ = i;                      /*   exp44[0] = alpha^0 */
    for(j = 0; j < 0x100; j++){
        a = (i>>4);                 /*   i *= alpha */
        b = (i&0xf);
        i = ((Mpy4(Mpy4(a,c),gf1)^Mpy4(a,d)^Mpy4(b,c))<<4)^Mpy4(Mpy4(a,c),gf0)^Mpy4(b,d);
        *pb++ = i;                  /*   exp44[j] = alpha^j */
        if(i == 0x01)break;}

    for(k = 8; k < 255; k += 1){    /* test mapping */
        j = exp08[k];               /*    j = GF(2^8) alpha^k */
        b44 = 0;                    /*    b44 = j mapped to GF(16^2) */
        for(i = 0; i < 8; i += 1){
            if(j&1){
                b44 ^= exp44[i];}
            j >>= 1;}
        if(exp44[k] != b44){        /*    if b44 != expected value */
            return(0);}}            /*      then not compatable */

    if((c == 1) && (d == 0)){       /* if GF(16^2) alpha == 1x + 0 */
        p44[p44idx].size = 3;       /*   save normal poly */
        p44[p44idx].data[0] = 1;
        p44[p44idx].data[1] = gf1;
        p44[p44idx].data[2] = gf0;
        p44idx += 1;}
    else{                           /* else */
        p4x[p4xidx].size = 3;       /*   save alternate poly */
        p4x[p4xidx].data[0] = 1;
        p4x[p4xidx].data[1] = gf1;
        p4x[p4xidx].data[2] = gf0;
        p4x[p4xidx].data[3] = c;
        p4x[p4xidx].data[4] = d;
        p4xidx += 1;}

    return(1);
}

/*----------------------------------------------------------------------*/
/*      PolyMpy(pPDst, pPSrc0, pPSrc1)  pPDst = pPSrc0 * pPSrc1         */
/*----------------------------------------------------------------------*/
static void PolyMpy(POLY *pPDst, POLY *pPSrc0, POLY *pPSrc1)
{
int i, j;

/*      init destination */

    pPDst->size = pPSrc0->size + pPSrc1->size - 1;
    memset(pPDst->data, 0, pPDst->size);

    for(j = 0; j < pPSrc0->size; j++){
        for(i = 0; i < pPSrc1->size; i++){
            pPDst->data[i+j] ^= Mpy4(pPSrc1->data[i], pPSrc0->data[j]);}}
}

/*----------------------------------------------------------------------*/
/*      InitTables                                                      */
/*----------------------------------------------------------------------*/
static void InitTables(void)            /* init tables */
{
BYTE *p0, *p1;
BYTE b0;

    b0 = 0x01;                           /* init exp4 table */
    p0 = exp4;
    for(p1 = p0+16; p0 < p1;){
        *p0++ = b0;
        b0 = Mp4x(b0, al4);}

    p0 = exp4;                          /* init log4 table */
    p1 = log4;
    *p1 = 0xf;
    for(b0 = 0; b0 < 15; b0 += 1)
        *(p1+*p0++) = b0;

    b0 = 0x01;                           /* init exp08 table */
    p0 = exp08;
    for(p1 = p0+256; p0 < p1;){
        *p0++ = b0;
        b0 = Mp8x(b0, al8);}
}

/*----------------------------------------------------------------------*/
/*      Mpy4    GF(16) multiply                                         */
/*----------------------------------------------------------------------*/
static BYTE Mpy4(BYTE m0, BYTE m1)      /* multiply */
{
int i;
    if(0 == m0 || 0 == m1)
        return(0);
    i = log4[m0] + log4[m1];
    if(i > 0xf)
        i -= 0xf;
    return(exp4[i]);
}

/*----------------------------------------------------------------------*/
/*      Div4    GF(16) divide                                           */
/*----------------------------------------------------------------------*/
static BYTE Div4(BYTE m0, BYTE m1)      /* divide */
{
int i;
    if(0 == m0)
        return(0);
    if(0 == m1){
        printf("Div4 - divide by zero\n");
        return(1);}
    i = log4[m0] - log4[m1];
    if(i < 0)
        i += 0xf;
    return(exp4[i]);
}

/*----------------------------------------------------------------------*/
/*      Mp4x(byt0, byt1)        returns product of byt0*byt1            */
/*                              no log or exp table used                */
/*----------------------------------------------------------------------*/
static BYTE Mp4x(BYTE byt0, BYTE byt1)
{
int i;
int product;
    product = 0;
    for(i = 0; i < 4; i++){
        product <<= 1;
        if(product & 0x10){
            product ^= gf4;}
        if(byt0 & 0x8u){
            product ^= byt1;}
        byt0 <<= 1;}
    return(product);
}

/*----------------------------------------------------------------------*/
/*      Mp8x(byt0, byt1)        returns product of byt0*byt1            */
/*                              no log or exp table used                */
/*----------------------------------------------------------------------*/
static BYTE Mp8x(BYTE byt0, BYTE byt1)
{
int i;
int product;
    product = 0;
    for(i = 0; i < 8; i++){
        product <<= 1;
        if(product & 0x100){
            product ^= gf8;}
        if(byt0 & 0x80u){
            product ^= byt1;}
        byt0 <<= 1;}
    return(product);
}
