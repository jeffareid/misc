/*----------------------------------------------------------------------*/
/*      crccor3.c   3 bit error corrrecting crc(1024,992)               */
/*                                                                      */
/*      Jeff Reid   2020JUN04  20:05                                    */
/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>

/*  1f1922815 =  787*557*465*3*3 */
/*      detects  7 err bits out of 1024 bits */ 
/*      corrects 3 err bits out of 1024 bits */
#define POLY 0xf1922815

typedef unsigned long long uint64_t;
typedef unsigned      long uint32_t;
typedef unsigned     short uint16_t;
typedef unsigned      char uint8_t;

/* comb(1024,1)+comb(1024,2)+comb(1024,3) = 178957824 */
#define C 178957824

/* N = number of bits (992 data + 32 crc) = 1024 */
#define N 1024

static uint32_t crctbb[N];      /* crc table, bit index */

static uint32_t crctbl[256];    /* table byte index */

static uint8_t  data[128];      /* test data */

static uint64_t *chktbl;        /* check table */

/*----------------------------------------------------------------------*/
/*      gentbl - generate crc tables                                    */
/*----------------------------------------------------------------------*/
void gentbl(void)
{
uint32_t crc;
uint32_t c;
int i;
    crc = 1;
    crctbb[N-1] = crc;                  /* generate crc bit  table */
    for(i = N-2; i >= 0; i--){
        crc = (crc<<1)^((0-(crc>>31))&POLY);
        crctbb[i] = crc;}
    for(c = 0; c < 0x100; c++){         /* generate crc byte table */
        crc = c<<24;
        for(i = 0; i < 8; i++){
            crc = (crc<<1)^((0-(crc>>31))&POLY);
        }
        crctbl[c] = crc;
    }
}

/*----------------------------------------------------------------------*/
/*      crc32c - generate crc                                           */
/*----------------------------------------------------------------------*/
uint32_t crc32c(uint8_t * bfr, size_t size)
{
uint32_t crc = 0;
    while(size--)
        crc = (crc << 8) ^ crctbl[(crc >> 24)^*bfr++];
    return(crc);
}

/*----------------------------------------------------------------------*/
/*      cmp64 - compare uint64_t    for qsort()                         */
/*----------------------------------------------------------------------*/
int cmp64(const void *j, const void *i)
{
    if((*(uint64_t *)j)<(*(uint64_t *)i))
        return -1;
    if((*(uint64_t *)j)>(*(uint64_t *)i))
        return  1;
    return 0;
}

/*----------------------------------------------------------------------*/
/*      InitCombination - init combination to first set - 1             */
/*----------------------------------------------------------------------*/
void InitCombination(int a[], int k, int n) {
    for(int i = 0; i < k; i++)
        a[i] = i;
    --a[k-1];     /* 1st call to NextCombination will return 1st set */
}

/*----------------------------------------------------------------------*/
/*      NextCombination - generate next combination                     */
/*----------------------------------------------------------------------*/
int NextCombination(int a[], int k, int n) {
int j = k - 1;
    while (j >= 0 && a[j] == n - k + j)
        --j;
    if (j == -1)
        return 0;
    ++a[j];
    for (int i = j + 1; i < k; ++i)
        a[i] = a[j] + i - j;
    return 1;
}

/*----------------------------------------------------------------------*/
/*      main                                                            */
/*----------------------------------------------------------------------*/
int main()
{
uint32_t crc;
uint32_t e[3];                          /* index to error bits */
int i, j, k;
    /* allocate 178,957,824*8 = 1,431,662,592 byte table */
    chktbl = (uint64_t *)malloc(C*sizeof(uint64_t));
    if(chktbl == (uint64_t *)0)
        goto exit0;
    gentbl();                           /* genreate crc tables */

    j = 0;
    InitCombination(e, 1, N);           /* add 1 bit error  to chktbl */
    while(NextCombination(e, 1, N))
        chktbl[j++] = (((uint64_t)(crctbb[e[0]]))<<32)
                            | (((uint64_t)e[0])<<20)|(((uint64_t)e[0])<<10)|(((uint64_t)e[0])<<0);
    InitCombination(e, 2, N);           /* add 2 bit errors to chktbl */
    while(NextCombination(e, 2, N))
        chktbl[j++] = (((uint64_t)(crctbb[e[0]]^crctbb[e[1]]))<<32)
                            | (((uint64_t)e[0])<<20)|(((uint64_t)e[1])<<10)|(((uint64_t)e[1])<<0);
    InitCombination(e, 3, N);           /* add 3 bit errors to chktbl */
    while(NextCombination(e, 3, N))
        chktbl[j++] = (((uint64_t)(crctbb[e[0]]^crctbb[e[1]]^crctbb[e[2]]))<<32)
                            | (((uint64_t)e[0])<<20)|(((uint64_t)e[1])<<10)|(((uint64_t)e[2])<<0);
    printf("chktbl created\n");
    /* sort chk table */
    qsort(chktbl, j, sizeof(uint64_t), cmp64);
    printf("chktbl sorted\n");
    for(i = 1; i < C; i++)              /* scan for duplicates */
        if((chktbl[i-1]&0xffffffff00000000ull) == (chktbl[i]&0xffffffff00000000ull))
            break;
    if(i == C)
        printf("no duplicates\n");
    else
        printf("   duplicates\n");
    
    e[0] = 320;                         /* simple test case */  
    e[1] = 640;
    e[2] = 960;
    data[e[0]>>3] |= 0x80u>>(e[0]&7);
    data[e[1]>>3] |= 0x80u>>(e[1]&7);
    data[e[2]>>3] |= 0x80u>>(e[2]&7);
/*  crc = crctbb[e[0]]^crctbb[e[1]]^crctbb[e[2]];  ** same crc as crc32c */
    crc = crc32c(data, 124);
    i = 0;
    j = C-1;
    do{                                 /* binary search chktbl for crc */
        k = (i+j)/2;
        if(((uint32_t)(chktbl[k]>>32) == crc))
            break;
        else if ((uint32_t)(chktbl[k]>>32) > crc)
            j = k-1;
        else
            i = k+1;
    }while(i <= j);
    if(i > j){
        printf("crc not found\n");
        goto exit0;}
    if(e[0] == (uint32_t)((chktbl[k]>>20)&0x3ff) &&
       e[1] == (uint32_t)((chktbl[k]>>10)&0x3ff) &&
       e[2] == (uint32_t)((chktbl[k]>> 0)&0x3ff))
       printf("passed\n");
    else
       printf("failed\n");
exit0:
    free(chktbl);
    return 0;
}
