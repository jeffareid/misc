/*----------------------------------------------------------------------*/
/*      ersp5c.c        RS(20,15) k=15 p=5 erasure code                 */
/*                                                                      */
/*      Copyright(c)    Jeff Reid   18OCT2025 18:15                     */
/*----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef unsigned char BYTE;
typedef unsigned int  WORD;
typedef unsigned long DWORD;
typedef unsigned long long QWORD;

static BYTE abExp[512];        /* antilog table */
static BYTE abLog[256];        /* log table */
static BYTE abInv[256];        /* inverse table */
static BYTE abTbl[4*64];       /* affine mpy table */
static BYTE *abDat;            /* allocated data */
static BYTE *pDat;             /* ptr, 64 byte bndry */
static BYTE abPoly[4] = {0xce, 0xe6, 0x00, 0x00};
/* poly = x^5 + ce x^4 + e6 x^3 + e6 x^2 + ce x + 01 */

static clock_t ctTimeStart;
static clock_t ctTimeStop;

extern void xenc(BYTE*, BYTE*, QWORD);

static BYTE gf_mul(BYTE b0, BYTE b1)
{
        if(b0 == 0 || b1 == 0)
            return(0);
        return(abExp[(int)abLog[b0]+(int)abLog[b1]]);
}

static BYTE gf_inv(BYTE b0)
{
        if(b0 == 0){
            printf("divide by zero\n");
            return(0);
        }
        return abInv[b0];
}

static BYTE gf_div(BYTE b0, BYTE b1)
{
        return gf_mul(b0, gf_inv(b1));
}

static BYTE gf_mul0(BYTE b0, BYTE b1)
{
int i;
int product;
        product = 0;
        for(i = 0; i < 8; i++){
            product <<= 1;
            if(product & 0x100u){
                product ^= 0x11du;}
            if(b0 & 0x80u){
                product ^= b1;}
            b0 <<= 1;}
        return((BYTE)product);
}

static void gf_init(void)
{
int i;
BYTE b, c, d, p;

    b = 1;                              /* init abExp[] */
    for(i = 0; i < 512; i++){
        abExp[i] = b;
        b = gf_mul0(b, 2);
    }

    abLog[0] = 0xffu;                   /* init abLog[] */
    for(i = 0; i < 255; i++)
        abLog[abExp[i]] = i;

    abInv[0] = 0xffu;                   /* init abInv[] */
    for(i = 1; i < 256; i++)
        abInv[i] = gf_div(1,i);

    for(i = 0; i < sizeof(abTbl);)      /* init abTbl */
    {
        d = abPoly[i/64];
        b = 0x80u;
        c = 7;
        do{
            p = gf_mul(b,d);
            abTbl[i+0] |= ((p>>7)&1)<<c;
            abTbl[i+1] |= ((p>>6)&1)<<c;
            abTbl[i+2] |= ((p>>5)&1)<<c;
            abTbl[i+3] |= ((p>>4)&1)<<c;
            abTbl[i+4] |= ((p>>3)&1)<<c;
            abTbl[i+5] |= ((p>>2)&1)<<c;
            abTbl[i+6] |= ((p>>1)&1)<<c;
            abTbl[i+7] |= ((p>>0)&1)<<c;
            c--;
        }while((b >>= 1) != 0);
        i += 8;
    }
}

#define NROW (20ull)
#define NCOL (2*1024*1024ull)

int main()
{
QWORD i;
    abDat = malloc(NROW*NCOL+64);
    pDat = (BYTE *)(((QWORD)abDat+63)&0xffffffffffffffc0ull);
    for (i = 0; i < NROW*NCOL; i++)
        pDat[i] = (BYTE)0;
#if 1
    for(i = 0; i < NCOL; i++)
        pDat[i] = 0x01u;
#endif
    gf_init();
    ctTimeStart = clock();
    for(i = 0; i < 1024; i++)
        xenc(abTbl, pDat, NCOL);
    ctTimeStop = clock();
    printf("# ticks: %d\n", ctTimeStop-ctTimeStart);
    free(abDat);
    return 0;
}