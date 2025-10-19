;-----------------------------------------------------------------------;
;       ersp5a.asm      RS(20,15) k=15 p=5 erasure code                 ;
;                                                                       ;
;       Copyright(c)    Jeff Reid   18OCT2025 18:15                     ;
;-----------------------------------------------------------------------;
        .data?
        align           16
        .data
        align           16
        PUBLIC  xenc
        .code
        align           16
;       xenc(rcx=abTbl, rdx=pDat, r8 = ncol)
;       rax, rcx, rdx, r8, r9, r10, r11 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5 free
;       k=15 p=5 poly = 01 x^5 + ce x^4 + e6 x^3 + e6 x^2 + ce x + 01
;                       x5       x4       x3       x2       x1     x0
xenc    proc
        mov             r9,r8                   ;rax = r8*19-64
        shl             r9,2
        lea             rax,[r9+r9*4]
        sub             rax,r8
        sub             rax,64
        mov             r9,r8                   ;r9 = ncol/64
        shr             r9,6
        vmovdqu64       zmm31,[rcx+  0]         ;mpy by ce
        vmovdqu64       zmm30,[rcx+ 64]         ;mpy by e6
        vmovdqu64       zmm29,zmm6              ;save zmm6

enc0:   vpxorq          zmm4,zmm4,zmm4          ;zero regs
        vpxorq          zmm3,zmm3,zmm3
        vpxorq          zmm2,zmm2,zmm2
        vpxorq          zmm1,zmm1,zmm1
        vpxorq          zmm0,zmm0,zmm0

;       logical values of zmm4 to zmm0 cycle during encode

        rept 3                                  ;encode 3*5 = 15 rows

        vpxorq          zmm4,zmm4,[rdx]         ;x0 ^= data
        vgf2p8affineqb  zmm5,zmm4,zmm31,0       ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm4,zmm30,0       ;z6 = x0*e6
        vpxorq          zmm3,zmm3,zmm5          ;x4 ^= z5
        vpxorq          zmm0,zmm0,zmm5          ;x1 ^= z5
        vpxorq          zmm2,zmm2,zmm6          ;x3 ^= z6
        vpxorq          zmm1,zmm1,zmm6          ;x2 ^= z6
        add             rdx,r8

        vpxorq          zmm3,zmm3,[rdx]         ;x0 ^= data
        vgf2p8affineqb  zmm5,zmm3,zmm31,0       ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm3,zmm30,0       ;z6 = x0*e6
        vpxorq          zmm2,zmm2,zmm5          ;x4 ^= z5
        vpxorq          zmm4,zmm4,zmm5          ;x1 ^= z5
        vpxorq          zmm1,zmm1,zmm6          ;x3 ^= z6
        vpxorq          zmm0,zmm0,zmm6          ;x2 ^= z6
        add             rdx,r8

        vpxorq          zmm2,zmm2,[rdx]         ;x0 ^= data
        vgf2p8affineqb  zmm5,zmm2,zmm31,0       ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm2,zmm30,0       ;z6 = x0*e6
        vpxorq          zmm1,zmm1,zmm5          ;x4 ^= z5
        vpxorq          zmm3,zmm3,zmm5          ;x1 ^= z5
        vpxorq          zmm0,zmm0,zmm6          ;x3 ^= z6
        vpxorq          zmm4,zmm4,zmm6          ;x2 ^= z6
        add             rdx,r8

        vpxorq          zmm1,zmm1,[rdx]         ;x0 ^= data
        vgf2p8affineqb  zmm5,zmm1,zmm31,0       ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm1,zmm30,0       ;z6 = x0*e6
        vpxorq          zmm0,zmm0,zmm5          ;x4 ^= z5
        vpxorq          zmm2,zmm2,zmm5          ;x1 ^= z5
        vpxorq          zmm4,zmm4,zmm6          ;x3 ^= z6
        vpxorq          zmm3,zmm3,zmm6          ;x2 ^= z6
        add             rdx,r8

        vpxorq          zmm0,zmm0,[rdx]         ;x0 ^= data
        vgf2p8affineqb  zmm5,zmm0,zmm31,0       ;z5 = x0*ce
        vgf2p8affineqb  zmm6,zmm0,zmm30,0       ;z6 = x0*e6
        vpxorq          zmm4,zmm4,zmm5          ;x4 ^= z5
        vpxorq          zmm1,zmm1,zmm5          ;x1 ^= z5
        vpxorq          zmm3,zmm3,zmm6          ;x3 ^= z6
        vpxorq          zmm2,zmm2,zmm6          ;x2 ^= z6
        add             rdx,r8

        endm

        vmovdqa64       [rdx],zmm4              ;par4 = x4
        add             rdx,r8
        vmovdqa64       [rdx],zmm3              ;par3 = x3
        add             rdx,r8
        vmovdqa64       [rdx],zmm2              ;par2 = x2
        add             rdx,r8
        vmovdqa64       [rdx],zmm1              ;par1 = x1
        add             rdx,r8
        vmovdqa64       [rdx],zmm0              ;par0 = x0

        sub             rdx,rax                 ;rdx = ptr to next 64 columns
        dec             r9                      ;loop till all columns done
        jnz             enc0

        vmovdqu64       zmm6,zmm29              ;restore zmm6
        ret
xenc    endp
        end
