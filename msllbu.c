/*  msllbu.c    merge sort linked list bottom up */
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

NODE * Merge(NODE *pNode0, NODE *pNode1)
{
NODE *pMrg;
NODE **ppMrg = &pMrg;
    if(pNode0 == NULL)                  /* if either list empty, */
        return pNode1;                  /*  return the other */
    if(pNode1 == NULL)
        return pNode0;
    while(1){                           /* merge the lists */
        if(pNode0->data <= pNode1->data){
            *ppMrg = pNode0;
            ppMrg = &(pNode0->next);
            pNode0 = *ppMrg;
            if(pNode0 == NULL){
                *ppMrg = pNode1;
                break;}}
        else{
            *ppMrg = pNode1;
            ppMrg = &(pNode1->next);
            pNode1 = *ppMrg;
            if(pNode1 == NULL){
                *ppMrg = pNode0;
                break;}}}
    return pMrg;
}


/* size of internal array */
#define ASZ 32
NODE *MergeSort(NODE * pNode)
{
NODE *apNode[ASZ] = {0};                /* array of lists */
NODE *pNext;
size_t i;
    if(pNode == NULL || pNode->next == NULL)
        return pNode;
    /* merge nodes into array */
    while(pNode){
        pNext = pNode->next;
        pNode->next = NULL;
        for(i = 0; i < ASZ && apNode[i] != NULL; i++){
            pNode = Merge(apNode[i], pNode);
            apNode[i] = NULL;}
        if(i == ASZ)                    /* don't go past end array */
            i--;
        apNode[i] = pNode;
        pNode = pNext;}
    /* merge array into a single list */
    pNode = NULL;
    for(i = 0; i < ASZ; i++)
        pNode = Merge(apNode[i], pNode);
    return pNode;
}

uint64_t rnd64()                        // random 64 bit integer
{
static uint64_t r = 0ull;
    r = r * 6364136223846793005ull + 1442695040888963407ull;
    return r;
}

NODE * TestSort(NODE *pNode)
{
NODE * pSort;
uint64_t r;
    pSort = pNode;
    while(pSort){
        r = rnd64();
        pSort->data =   r;
        pSort = pSort->next;
#if 0
        pSort->data = ++r;
        pSort = pSort->next;
        pSort->data = ++r;
        pSort = pSort->next;
        pSort->data = ++r;
        pSort = pSort->next;
        pSort->data = ++r;
        pSort = pSort->next;
        pSort->data = ++r;
        pSort = pSort->next;
        pSort->data = ++r;
        pSort = pSort->next;
        pSort->data = ++r;
        pSort = pSort->next;
#endif
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

    pSort = MergeSort(pNode);

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
#define COUNT (16*1024*1024)
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

    pNode = TestSort(pNode);
    pNode = TestSort(pNode);

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
