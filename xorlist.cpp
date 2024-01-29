//----------------------------------------------------------------------//
//      xorlist.cpp     xor linked list                                 //
//----------------------------------------------------------------------//
#include <ctime>
#include <iostream>
#include <math.h>

typedef unsigned long long uint64_t;

class XORNODE{
public:
    size_t   link;
    uint64_t data;
    XORNODE(){
        this->link = 0;
        this->data = 0;
    }
    XORNODE(uint64_t d){
        this->link = 0;
        this->data = d;
    }
    ~XORNODE(){
    }
};

class XORITER{
public:
    XORNODE *prv;
    XORNODE *cur;
    XORITER(){
        this->prv = NULL;
        this->cur = NULL;
    }
    ~XORITER(){
    }
    XORITER next(){
        XORNODE *prv = this->prv;
        XORNODE *cur = this->cur;
        this->prv = cur;
        this->cur = ((XORNODE *)((size_t)(prv)^(size_t)(cur->link)));
        return *this;
    }
    XORITER prev(){
        XORNODE *prv = this->prv;
        XORNODE *cur = this->cur;
        this->prv = ((XORNODE *)((size_t)(cur)^(size_t)(prv->link)));
        this->cur = prv;
        return *this;
    }
};

class XORLIST{
    XORNODE *bbb;
    XORNODE *eee;
    XORNODE dmy;
    size_t size;
public:
    XORLIST(){
        bbb = NULL;
        eee = &dmy;
        size = 0;
    }
    ~XORLIST(){
    }

    XORITER begin(){
        XORITER itr;
        itr.prv = this->eee;
        itr.cur = this->bbb;
        return itr;
    }

    XORITER end(){
        XORITER itr;
        XORNODE *lst = (XORNODE *)((size_t)(eee->link)^(size_t)(bbb));
        itr.prv = lst;
        itr.cur = this->eee;
        return itr;
    }

    void push_back(uint64_t d){
        size++;
        XORNODE* pnw = new XORNODE;
        pnw->data = d;
        if(bbb == NULL){
            bbb = pnw;
            return;
        }
        XORNODE *lst = (XORNODE *)((size_t)(eee->link)^(size_t)(bbb));
        lst->link ^= ((size_t)(eee) ^ (size_t)(pnw));
        pnw->link  = ((size_t)(lst) ^ (size_t)(eee));
        eee->link  = ((size_t)(pnw) ^ (size_t)(bbb));
    }

    XORITER insert(XORITER pos, uint64_t d){
        size++;
        XORNODE* pnw = new XORNODE;
        pnw->data = d;
        if(bbb == NULL){                // ignore pos if empty list
            bbb = pnw;
            pos.prv = eee;
            pos.cur = pnw;
            return pos;
        }
        XORNODE *prv = pos.prv;
        XORNODE *cur = pos.cur;
        if(prv == eee)
            bbb = pnw;
        prv->link ^= ((size_t)(cur)^(size_t)(pnw));
        pnw->link  = ((size_t)(prv)^(size_t)(cur));
        cur->link ^= ((size_t)(prv)^(size_t)(pnw));
        pos.cur = pnw;
        return pos;
    }

    XORITER erase(XORITER pos){
        if(bbb == NULL)                 // ignore if empty list
            return pos;
        size--;
        XORNODE *prv = pos.prv;
        XORNODE *cur = pos.cur;
        XORNODE *nxt = ((XORNODE *)((size_t)(prv)^(size_t)(cur->link)));
        prv->link ^= ((size_t)(cur)^(size_t)(nxt));
        nxt->link ^= ((size_t)(cur)^(size_t)(prv));
        delete cur;
        pos.cur = nxt;
        return pos;
    }

    void reverse(){
        bbb = (XORNODE *)((size_t)(eee->link)^(size_t)(bbb));
    }

    void splice(XORITER pos, XORLIST &othr){
        if(othr.bbb == NULL)
            return;
        size += othr.size;
        XORNODE *prv = pos.prv;
        XORNODE *cur = pos.cur;
        XORNODE *ols = (XORNODE *)((size_t)(othr.eee->link)^(size_t)(othr.bbb));
        prv->link      ^= ((size_t)(cur)^(size_t)(othr.bbb));
        cur->link      ^= ((size_t)(prv)^(size_t)(ols));
        othr.bbb->link ^= ((size_t)(prv)^(size_t)(othr.eee));
        ols->link      ^= ((size_t)(cur)^(size_t)(othr.eee));
        othr.bbb = NULL;
        othr.eee->link = 0;
    }

    void splice(XORITER &lft, XORITER &rgt){
        XORNODE *lfp = lft.prv;
        XORNODE *lfc = lft.cur;
        XORNODE *lfn = (XORNODE *)((size_t)(lfp)^(size_t)(lfc->link));
        XORNODE *rgp = rgt.prv;
        XORNODE *rgc = rgt.cur;
        XORNODE *rgn = (XORNODE *)((size_t)(rgp)^(size_t)(rgc->link));
        if(lfp == eee)
            bbb = rgc;
        lfp->link ^= ((size_t)(lfc)^(size_t)(rgc));
        rgc->link  = ((size_t)(lfp)^(size_t)(lfc));
        lfc->link  = ((size_t)(rgc)^(size_t)(lfn));
        rgp->link ^= ((size_t)(rgc)^(size_t)(rgn));
        rgn->link ^= ((size_t)(rgc)^(size_t)(rgp));
        lft.prv = rgc;
        rgt.cur = rgn;
        if(lfc == rgp)
            rgt.prv = lfc;
    }

    XORITER merge(XORITER li, XORITER ri, XORITER &ei){
        XORITER ni = li;
        if(ri.cur->data < li.cur->data)
            ni.cur = ri.cur;
        while(1){
            if(ri.cur->data < li.cur->data){
                splice(li, ri);
                if (ri.cur == ei.cur) {
                    ei.prv = ri.prv;
                    return ni;
                }
            } else {
                li.next();
                if(li.cur == ri.cur)
                    return ni;
            }
        }
    }

    XORITER sortr(XORITER &li, size_t sz)
    {
        XORITER ri;
        XORITER ei;
        if (sz == 1){
            ei = li;
            return ei.next();
        }
        ri = sortr(li, sz-sz/2);
        ei = sortr(ri, sz/2);
        li = merge(li, ri, ei);
        return ei;
    }

    void sort()
    {
        if(size < 2)                    // return if nothing to do
            return;
        sortr(begin(), size);
    }

    void show(){
        XORITER itr = this->begin();
        XORITER stp = this->end();
        if(itr.cur == NULL)
            return;
        while(itr.cur != stp.cur){
            std::cout << itr.cur->data << ' ';
            itr.next();
        }
        std::cout << std::endl;
#if 0
        itr = this->end();
        stp = this->begin();
        while(itr->cur != stp->cur){
            itr->prev();
            std::cout << itr->cur->data << ' ';
        }
        std::cout << std::endl;
#endif
    }
};

#if 1

uint64_t rnd64()                        // random 64 bit integer
{
static uint64_t r = 0ull;
    r = r * 6364136223846793005ull + 1442695040888963407ull;
    return r;
}

#define COUNT (4*1024*1024)

int main(void)
{
XORLIST list;
XORITER iter;
XORITER iend;
uint64_t r;
size_t i;
    for(i = 0; i < COUNT; i++)
        list.push_back(rnd64());
    list.sort();
    iter = list.begin();
    iend = list.end();
    r = iter.cur->data;
    iter.next();
    while(iter.cur != iend.cur){
        if(r > iter.cur->data)
            break;
        r = iter.cur->data;
        iter.next();
    }
    if(iter.cur != iend.cur)
        std::cout << "failed" << std::endl;
    else
        std::cout << "passed" << std::endl;
    return(0);
}

#else

int main(void)
{
XORLIST list;
XORLIST othr;
XORITER iter;
XORITER lft, rgt;
    list.push_back(1);
    list.push_back(2);
    list.push_back(5);
    iter = list.begin();
    iter.next();
    iter.next();
    iter = list.insert(iter, 4);
    iter = list.insert(iter, 3);
    list.show();
    othr.push_back(9);
    othr.push_back(8);
    othr.push_back(7);
    othr.push_back(6);
    othr.show();
    othr.reverse();
    list.splice(list.end(), othr);
    list.show();
    list.reverse();
    list.show();
    list.sort();
    list.show();
    return(0);
}

#endif