/*  msllrc.c    merge sort linked list recursive */
#include <stdio.h>
#include <stdlib.h>

typedef unsigned long long uint64_t;

/* if QP != 0, use queryperformance for timer */
#define QP 1

#if QP
#include <math.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
typedef LARGE_INTEGER LI64;
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
static clock_t ctTimeStart;     /* clock values */
static clock_t ctTimeStop;
#endif


typedef struct NODE_{                       /* base node struct */
    struct NODE_ * next;
    uint64_t       data;
}NODE;

typedef struct{
    NODE *first;
    NODE *last;
}NODEPAIR;

NODEPAIR Merge(NODE *prev, NODE *p0, int n0,
                   NODE *p1, int n1)
{
NODEPAIR np;
NODE *first;
int i, j;
    first = p0->data <= p1->data ? p0 : p1;
    i = j = 0;
    while(1){
        if(p0->data <= p1->data){       /* if p0 < p1 */
            prev->next = p0;            /*  move p0 */
            prev = p0;
            p0 = p0->next;
            if(++i < n0)                /*  if not end run 0 */
                continue;               /*   continue back to while */
            prev->next = p1;            /*  else link run 1 */
            while(++j < n1)
                p1 = p1->next;
            prev = p1;
            p1 = p1->next;
            break;
        } else {
            prev->next = p1;            /*  move p1 */
            prev = p1;
            p1 = p1->next;
            if(++j < n1)                /*  if not end run 1 */
                continue;               /*   continue back to while */
            prev->next = p0;            /*  else link run 0 */
            while(++i < n0)
                p0 = p0->next;
            prev = p0;
            p0 = p0->next;
            break;
        }
    }
    np.first = first;
    np.last  = prev;
    prev->next = p1;
    return np;
}

NODEPAIR MergeSortR(NODE *prev, NODE *p0, int n0)
{
NODEPAIR np;
NODE *p1;
int n1;
    if (n0 <= 1){
        np.first = np.last = p0;
        return np;
    }
    n1 = n0;
    n0 >>= 1;
    n1 -= n0;
    np = MergeSortR(prev, p0, n0);
    p0 = np.first;
    p1 = np.last->next;
    np = MergeSortR(np.last, p1, n1);
    p1 = np.first;
    np = Merge(prev, p0, n0, p1, n1);
    return np;
}

NODE * MergeSort(NODE *p0, int n0)
{
NODE prevnode = {0,0};
NODE * prev = &prevnode;
    return MergeSortR(prev, p0, n0).first;
}

uint64_t rnd64()                        // random 64 bit integer
{
static uint64_t r = 0ull;
    r = r * 6364136223846793005ull + 1442695040888963407ull;
    return r;
}

NODE * TestSort(NODE *pNode, int n)
{
NODE * pSort;
    pSort = pNode;
    while(pSort){
        pSort->data = rnd64();
        pSort = pSort->next;
    }

#if QP
    QueryPerformanceFrequency(&liQPFrequency);
    dQPFrequency = (double)liQPFrequency.QuadPart;
    timeBeginPeriod(250);                   /* set ticker to 4 hz */
    Sleep(250);                             /* wait for it to settle */
    QueryPerformanceCounter(&liStartTime);
#else
    ctTimeStart = clock();
#endif

    pSort = MergeSort(pNode, n);

#if QP
    QueryPerformanceCounter(&liStopTime);
    dStartTime = (double)liStartTime.QuadPart;
    dStopTime  = (double)liStopTime.QuadPart;
    dElapsedTime = (dStopTime - dStartTime) / dQPFrequency;
    timeEndPeriod(250);                     /* restore ticker to default */
    printf("# of seconds %f\n", dElapsedTime);
#else
    ctTimeStop = clock();
    printf("# of ticks %u\n", ctTimeStop - ctTimeStart);
#endif

    return pSort;
}

/* # of nodes */
#define COUNT (16*1024*1024+1)
int main()                                  /* main */
{
NODE * pMem = (NODE *)malloc(COUNT*sizeof(NODE));
NODE * pNode;
uint64_t r;
size_t i;
    pNode = pMem;
    for(i = 0; i < COUNT; i++)
        pNode[i].next = pNode+i+1;
    pNode[i-1].next = NULL;

    pNode = TestSort(pNode, COUNT);
    pNode = TestSort(pNode, COUNT);

    r = pNode->data;
    for(i = 1; i < COUNT; i++){
        pNode = pNode->next;
        if(r > pNode->data)
            break;
        r = pNode->data;
    }
    if(i == COUNT)
        printf("passed\n");
    else
        printf("failed\n");

    free(pMem);
    return(0);
}
