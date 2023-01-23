/*----------------------------------------------------------------------*/
/*      msbu4c32.c      hybrid insertion + merge sort + cmp func        */
/*                                                                      */
/*      Jeff Reid       2022AUG12 09:15                                 */
/*      modified        2023JAN23 12:50                                 */
/*----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

typedef int CMPFUNC (const void *a, const void *b);
#if 1
#define cmp(a,b) (*(a) > *(b))
#endif

int cmp_int(const void * a, const void * b)
{
    return *(int *)a > *(int *)b;
}

#define QP 0        /* if != 0, use queryperformance for timer */

#if QP
#include <math.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
typedef LARGE_INTEGER LI64;
#else
#include <time.h>
#endif

#if QP
static LI64     liQPFrequency;  /* cpu counter values */
static LI64     liStartTime;
static LI64     liStopTime;
static double   dQPFrequency;
static double   dStartTime;
static double   dStopTime;
static double   dElapsedTime;
#else
clock_t ctTimeStart;            /* clock values */
clock_t ctTimeStop;
#endif

size_t GetPassCount(size_t n);

void MergeSort(int a[], size_t n, CMPFUNC cmp)
{
int *p0r;                                   /* ptr to current element run 0 */
int *p0e;                                   /* ptr to end             run 0 */
int *p1r;                                   /* ptr to current element run 1 */
int *p1e;                                   /* ptr to end             run 1 */
int *p2r;                                   /* ptr to current element run 2 */
int *p2e;                                   /* ptr to end             run 2 */
int *p3r;                                   /* ptr to current element run 3 */
int *p3e;                                   /* ptr to end             run 3 */
int *pax;                                   /* ptr to a[] or b[] */
int *pbx;                                   /* ptr to b[] or a[] */
int * b;                                    /* second array */
size_t rsz;                                 /* run size  */
    if(n < 4){                              /* check if size < 4 */
        int t;
        if(n < 2)
            return;
        if(n == 2){
            if(cmp(&a[0], &a[1])){
                t    = a[0];
                a[0] = a[1];
                a[1] = t;
            }
            return;
        }
        if(n == 3){
            if(cmp(&a[0], &a[2])){
                t    = a[0];
                a[0] = a[2];
                a[2] = t;
            }
            if(cmp(&a[0], &a[1])){
                t    = a[0];
                a[0] = a[1];
                a[1] = t;
            }
            if(cmp(&a[1], &a[2])){
                t    = a[1];
                a[1] = a[2];
                a[2] = t;
            }
            return;
        }
    }
    b = malloc(n * sizeof(int));

    /* set run size so merge sort is even number of passes */
    rsz = ((GetPassCount(n) & 1) != 0) ? 64 : 16;
    {                                       /* insertion sort */
        size_t l, r;
        size_t i, j;
        int t;
        for (l = 0; l < n; l = r) {
            r = l + rsz;
            if (r > n)r = n;
            l--;
            for (j = l + 2; j < r; j++) {
                t = a[j];
                i = j-1;
                while(i != l && cmp(&a[i], &t)){
                    a[i+1] = a[i];
                    i--;
                }
                a[i+1] = t;
            }
        }
    }

    while(rsz < n){                         /* merge sort */
        pbx = &b[0];
        pax = &a[0];
        while(pax < &a[n]){
            p0e = rsz + (p0r = pax);
            if(p0e >= &a[n]){
                p0e = &a[n];
                goto cpy10;}
            p1e = rsz + (p1r = p0e);
            if(p1e >= &a[n]){
                p1e = &a[n];
                goto mrg201;}
            p2e = rsz + (p2r = p1e);
            if(p2e >= &a[n]){
                p2e = &a[n];
                goto mrg3012;}
            p3e = rsz + (p3r = p2e);
            if(p3e >= &a[n])
                p3e = &a[n];
            /* 4 way merge */
            while(1){
                if(cmp(p1r, p0r)){
                    if(cmp(p3r, p2r)){
                        if(cmp(p2r, p0r)){
mrg40:                      *pbx++ = *p0r++;    /* run 0 smallest */
                            if(p0r < p0e)       /* if not end run continue */
                                continue;
                            goto mrg3123;       /* merge 1,2,3 */
                        } else {
mrg42:                      *pbx++ = *p2r++;    /* run 2 smallest */
                            if(p2r < p2e)       /* if not end run continue */
                                continue;
                            goto mrg3013;       /* merge 0,1,3 */
                        }
                    } else {
                        if(cmp(p3r, p0r)){
                            goto mrg40;         /* run 0 smallext */
                        } else {
mrg43:                      *pbx++ = *p3r++;    /* run 3 smallest */
                            if(p3r < p3e)       /* if not end run continue */
                                continue;
                            goto mrg3012;       /* merge 0,1,2 */
                        }
                    }
                } else {
                    if(cmp(p3r, p2r)){
                        if(cmp(p2r, p1r)){
mrg41:                      *pbx++ = *p1r++;    /* run 1 smallest */
                            if(p1r < p1e)       /* if not end run continue */
                                continue;
                            goto mrg3023;       /* merge 0,2,3 */
                        } else {
                            goto mrg42;         /* run 2 smallest */
                        }
                    } else {
                        if(cmp(p3r, p1r)){
                            goto mrg41;         /* run 1 smallest */
                        } else {
                            goto mrg43;         /* run 3 smallest */
                        }
                    }
                }
            }
            /* 3 way merge */
mrg3123:    p0r = p1r;
            p0e = p1e;
mrg3023:    p1r = p2r;
            p1e = p2e;
mrg3013:    p2r = p3r;
            p2e = p3e;
mrg3012:    while(1){
                if(cmp(p1r, p0r)){
                    if(cmp(p2r, p0r)){
                        *pbx++ = *p0r++;        /* run 0 smallest */
                        if(p0r < p0e)           /* if not end run continue */
                            continue;
                        goto mrg212;            /* merge 1,2 */
                    } else {
mrg32:                  *pbx++ = *p2r++;        /* run 2 smallest */
                        if(p2r < p2e)           /* if not end run continue */
                            continue;
                        goto mrg201;            /* merge 0,1 */
                    }
                } else {
                    if(cmp(p2r, p1r)){
                        *pbx++ = *p1r++;        /* run 1 smallest */
                        if(p1r < p1e)           /* if not end run continue */
                            continue;
                        goto mrg202;            /* merge 0,2 */
                    } else {
                        goto mrg32;             /* run 2 smallest */
                    }
                }
            }
mrg212:     p0r = p1r;
            p0e = p1e;
mrg202:     p1r = p2r;
            p1e = p2e;
            /* 2 way merge */
mrg201:     while(1){
                if(cmp(p1r, p0r)){
                    *pbx++ = *p0r++;            /* run 0 smallest */
                    if(p0r < p0e)               /* if not end run continue */
                        continue;
                    goto cpy11;
                } else {
                    *pbx++ = *p1r++;            /* run 1 smallest */
                    if(p1r < p1e)               /* if not end run continue */
                        continue;
                    goto cpy10;
                }
            }
cpy11:      p0r = p1r;
            p0e = p1e;
            /* 1 way copy */
cpy10:      while (1) {
                *pbx++ = *p0r++;                /* copy element */
                if (p0r < p0e)                  /* if not end of run continue */
                    continue;
                break;
            }
            pax += rsz << 2;                    /* setup for next set of runs */
        }
        pax = a;                                /* swap ptrs */
        a = b;
        b = pax;
        rsz <<= 2;                              /* quadruple run size */
    }
    free(b);
}

size_t GetPassCount(size_t n)                   /* return # passes */
{
    size_t i = 0;
    size_t s;
    for(s = 1; s < n; s <<= 2)
        i += 1;
    return(i);
}

int rnd32()
{
static unsigned int r = 0;
    r = r*1664525 + 1013904223;
    return (int)r;
}

#define COUNT (0x1000000)

int main()
{
int * a = malloc(COUNT * sizeof(int));
size_t i;
    for(i = 0; i < COUNT; i++)
        a[i] = rnd32();
#if QP
    QueryPerformanceFrequency(&liQPFrequency);
    dQPFrequency = (double)liQPFrequency.QuadPart;
    timeBeginPeriod(256);                   /* set ticker to 4 hz */
    Sleep(512);                             /* wait for it to settle */
    QueryPerformanceCounter(&liStartTime);
#else
    ctTimeStart = clock();
#endif
    MergeSort(a, COUNT, cmp_int);
#if QP
    QueryPerformanceCounter(&liStopTime);
    dStartTime = (double)liStartTime.QuadPart;
    dStopTime  = (double)liStopTime.QuadPart;
    dElapsedTime = (dStopTime - dStartTime) / dQPFrequency;
    timeEndPeriod(256);                     /* restore ticker to default */
    printf("# of seconds %f\n", dElapsedTime);
#else
    ctTimeStop = clock();
    printf("# of ticks %d\n", ctTimeStop - ctTimeStart);
#endif
    for(i = 1; i < COUNT; i++){
        if(a[i] < a[i-1])
            break;
    }
    if(i == COUNT)
        printf("passed\n");
    else
        printf("failed\n");
    free(a);
    return 0;
}
