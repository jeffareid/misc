;-----------------------------------------------------------------------;
;       ersp5a.asm      RS(20,15) k=15 p=5 erasure code                 ;
;                                                                       ;
;       Copyright(c)    Jeff Reid   18OCT2025 18:30                     ;
;-----------------------------------------------------------------------;
        .data?
        align           16
        .data
        align           16
        PUBLIC  xenc
        .code
        align           16
;       xenc(rcx=abTbl, rdx=pDat, r8 = ncol)
;       k=15 p=5 poly = 01 x^5 + ce x^4 + e6 x^3 + e6 x^2 + ce x + 01
;                                x4       x3       x2       x1     x0
;       rax, rcx, rdx, r8, r9, r10, r11 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5 free
xenc    proc
        vmovdqa64       zmm29,zmm9              ;save
        vmovdqa64       zmm28,zmm8
        vmovdqa64       zmm27,zmm7
        vmovdqa64       zmm26,zmm6
        vmovdqu64       zmm9,[rcx+  0]          ;mpy by ce
        vmovdqu64       zmm8,[rcx+ 64]          ;mpy by e6
        mov             rax,r8                  ;rax = r8*15-64
        shl             rax,4
        sub             rax,r8
        sub             rax,64
        mov             r9,r8                   ;r9 = ncol/64
        shr             r9,6
        lea             r10,[r8+r8*2]           ;r10 = r8*3
        lea             r11,[r8+r8*4]           ;r11 = r8*5

enc0:   vmovdqa64       zmm7,[rdx]              ;preload data
        vpxorq          zmm3,zmm3,zmm3          ;z3 = 0
        vpxorq          zmm2,zmm2,zmm2          ;z2 = 0
        vpxorq          zmm1,zmm1,zmm1          ;z1 = 0
        vmovdqa64       zmm4,zmm7               ;z4 = data
        vpxorq          zmm0,zmm0,zmm0          ;z0 = 0

;       logical values of zmm4 to zmm0 cycle during encode

encd    =               3
        rept            3                       ;encode 3*5 = 15 rows
encd    =               encd-1
        vmovdqa64       zmm7,[rdx+r8]           ;preload data
        vgf2p8affineqb  zmm5,zmm4,zmm9,0        ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm4,zmm8,0        ;z6 = x0*e6
        vpxorq          zmm0,zmm0,zmm5          ;x1 ^= z5
        vpternlogq      zmm3,zmm5,zmm7,096h     ;x4 ^= z5^data
        vpxorq          zmm1,zmm1,zmm6          ;x2 ^= z6
        vpxorq          zmm2,zmm2,zmm6          ;x3 ^= z6

        vmovdqa64       zmm7,[rdx+r8*2]         ;preload data
        vgf2p8affineqb  zmm5,zmm3,zmm9,0        ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm3,zmm8,0        ;z6 = x0*e6
        vpxorq          zmm4,zmm4,zmm5          ;x1 ^= z5
        vpternlogq      zmm2,zmm5,zmm7,096h     ;x4 ^= z5^data
        vpxorq          zmm0,zmm0,zmm6          ;x2 ^= z6
        vpxorq          zmm1,zmm1,zmm6          ;x3 ^= z6

        vmovdqa64       zmm7,[rdx+r10]          ;preload data
        vgf2p8affineqb  zmm5,zmm2,zmm9,0        ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm2,zmm8,0        ;z6 = x0*e6
        vpxorq          zmm3,zmm3,zmm5          ;x1 ^= z5
        vpternlogq      zmm1,zmm5,zmm7,096h     ;x4 ^= z5^data
        vpxorq          zmm4,zmm4,zmm6          ;x2 ^= z6
        vpxorq          zmm0,zmm0,zmm6          ;x3 ^= z6

        vmovdqa64       zmm7,[rdx+r8*4]         ;preload data
        vgf2p8affineqb  zmm5,zmm1,zmm9,0        ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm1,zmm8,0        ;z6 = x0*e6
        vpxorq          zmm2,zmm2,zmm5          ;x1 ^= z5
        vpternlogq      zmm0,zmm5,zmm7,096h     ;x4 ^= z5^data
        vpxorq          zmm3,zmm3,zmm6          ;x2 ^= z6
        vpxorq          zmm4,zmm4,zmm6          ;x3 ^= z6

        if              encd
        vmovdqa64       zmm7,[rdx+r11]          ;preload data
        endif
        vgf2p8affineqb  zmm5,zmm0,zmm9,0        ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm0,zmm8,0        ;z6 = x0*e6
        vpxorq          zmm1,zmm1,zmm5          ;x1 ^= z5
        if              encd
        vpternlogq      zmm4,zmm5,zmm7,096h     ;x4 ^= z5^data
        else
        vpxorq          zmm4,zmm4,zmm5          ;x4 ^= z5
        endif
        vpxorq          zmm2,zmm2,zmm6          ;x2 ^= z6
        vpxorq          zmm3,zmm3,zmm6          ;x3 ^= z6
        add             rdx,r11                 ;rdx += ncol*5

        endm

        vmovdqa64       [rdx],zmm4              ;par4 = x4
        vmovdqa64       [rdx+r8],zmm3           ;par3 = x3
        vmovdqa64       [rdx+r8*2],zmm2         ;par2 = x2
        vmovdqa64       [rdx+r10],zmm1          ;par1 = x1
        vmovdqa64       [rdx+r8*4],zmm0         ;par0 = x0

        sub             rdx,rax                 ;rdx = ptr to next 64 columns
        dec             r9                      ;loop till all columns done
        jnz             enc0

        vmovdqa64       zmm6,zmm26              ;restore
        vmovdqa64       zmm7,zmm27
        vmovdqa64       zmm8,zmm28
        vmovdqa64       zmm9,zmm29
        ret
xenc    endp
        end
