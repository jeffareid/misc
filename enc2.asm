        page    80,128
        title   enc2
;-------------------------------------------------------;
;       enc2    286 enc code                            ;
;                                                       ;
;       Copyright(c)	Jeff Reid    01MAR88 11:00      ;
;-------------------------------------------------------;
;       equates                                         ;
;-------------------------------------------------------;
        .286
        .sall
        .model  small,c
;-------------------------------------------------------;
;       unitialized data section                        ;
;-------------------------------------------------------;
        .data?
log2    db      256 dup (?)     ;log2 table
alog2   db      256 dup (?)     ;alog2 table
mxc0    db      256 dup (?)     ;mpy by c0 table
bfr     db      32768 dup (?)   ;bfr
;-------------------------------------------------------;
;       stack                                           ;
;-------------------------------------------------------;
        .stack  1024
;-------------------------------------------------------;
;       data section                                    ;
;-------------------------------------------------------;
        .data
;-------------------------------------------------------;
;       code section                                    ;
;-------------------------------------------------------;
        .code
        assume  cs:@code,ds:@data,es:nothing,ss:nothing

;-------------------------------------------------------;
;       enc     encode data                             ;
;               in: 2(sp) = virt adr bfr                ;
;-------------------------------------------------------;
r0      equ     al              ;define remainder bytes
r1      equ     cl
r2      equ     dl
;
enc     proc    near
        push    bp              ;set up
        mov     bp,sp
        push    es
        push    di
        xor     bx,bx
        mov     si,1024
        les     di,4[bp]        ;es:di = adr
;
enc0:   xor     ax,ax
        xor     cx,cx
        xor     dx,dx
encd    =       0
        xor     r1,es:[di+encd*1024]    ;do 1st 2 bytes
        mov     bl,r1
        xor     r0,mxc0[bx]
        xor     r2,mxc0[bx]
encd    =       encd+1
        xor     r0,es:[di+encd*1024]
        mov     bl,r0
        xor     r2,mxc0[bx]
        xor     r1,mxc0[bx]
encd    =       encd+1

        rept    9                       ;do remaining 27 bytes
        xor     r2,es:[di+encd*1024]
        mov     bl,r2
        xor     r1,mxc0[bx]
        xor     r0,mxc0[bx]
encd    =       encd+1
        xor     r1,es:[di+encd*1024]
        mov     bl,r1
        xor     r0,mxc0[bx]
        xor     r2,mxc0[bx]
encd    =       encd+1
        xor     r0,es:[di+encd*1024]
        mov     bl,r0
        xor     r2,mxc0[bx]
        xor     r1,mxc0[bx]
encd    =       encd+1
        endm

        mov     es:[di+29*1024],r2      ;store encoded data
        mov     es:[di+30*1024],r1
        mov     es:[di+31*1024],r0
        inc     di                      ;adv adr
        dec     si
        jnz     enc0
        pop     di
        pop     es
        pop     bp
        ret
enc     endp
;-------------------------------------------------------;
;       tbli    initialize table                        ;
;-------------------------------------------------------;
tbli    proc    near
        mov     bx,offset alog2 ;init alog2 table
        mov     al,01
tbli0:  mov     [bx],al
        inc     bx
        add     al,al
        jnc     tbli1
        xor     al,87h
tbli1:  cmp     al,01
        jne     tbli0
;
        xor     bx,bx           ;init log2 table
        mov     si,offset alog2
        mov     log2,0ffh
        xor     ax,ax
tbli2:  mov     bl,[si]
        mov     log2[bx],al
        inc     si
        inc     ax
        cmp     ax,255
        jne     tbli2
;
        mov     si,offset mxc0  ;init mpy by c0 table
        mov     cx,0c000h
tbli3:  call    mpy
        mov     [si],al
        inc     si
        inc     cl
        jne     tbli3
        ret
tbli    endp
;-------------------------------------------------------;
;       mpy     ax = ch*cl                              ;
;-------------------------------------------------------;
mpy     proc    near
        xor     ax,ax           ;return 0 if 0
        or      cl,cl
        je      mpy1
        or      ch,ch
        je      mpy1
        xor     bx,bx           ;ax=log[ch]+log[cl]
        mov     bl,ch
        mov     al,log2[bx]
        mov     bl,cl
        mov     bl,log2[bx]
        add     ax,bx
        cmp     ax,255          ;ax = mod(ax, 255)
        jb      mpy0
        sub     ax,255
mpy0:   xchg    ax,bx           ;ax = alog2[ax]
        xor     ax,ax
        mov     al,alog2[bx]
mpy1:   ret
mpy     endp

;-------------------------------------------------------;
;       main                                            ;
;-------------------------------------------------------;
main    proc    far
        mov     ax,@data
        mov     ds,ax
        mov     es,ax
        call    tbli
        lea     di,bfr
        mov     cx,32768
        xor     ax,ax
        rep     stosb
        lea     di,bfr
        mov     cx,1024
        mov     al,01h
        rep     stosb
        lea     di,bfr
        push    es
        push    di
        call    enc
        add     sp,4
        mov     ax,4c00h
        int     21h
main    endp

        end     main

