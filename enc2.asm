        page    80,128
        title   enc2
;-------------------------------------------------------;
;       enc2    286 enc code                            ;
;                                                       ;
;       Jeff Reid       01MAR88 11:00                   ;
;-------------------------------------------------------;
;       equates                                         ;
;-------------------------------------------------------;
        .286
        .sall
        .model  small,c
stksz   equ     00400h          ;initial stack size
wrksz   equ     00800h          ;work | stack size
;-------------------------------------------------------;
;       unitialized data section                        ;
;-------------------------------------------------------;
        .data?
;-------------------------------------------------------;
;       stack                                           ;
;-------------------------------------------------------;
        .stack  stksz           ;initial stack size
;-------------------------------------------------------;
;       data section                                    ;
;-------------------------------------------------------;
        .data
;-------------------------------------------------------;
;       work | stack segment                            ;
;-------------------------------------------------------;
exp2    equ     00000h          ;exp2 table
log2    equ     00200h          ;log2 table
mxc0    equ     00300h          ;mpy by c0 table
;-------------------------------------------------------;
;       bfr segment                                     ;
;-------------------------------------------------------;
bfr     equ     00000h          ;buffer
;-------------------------------------------------------;
;       code section                                    ;
;-------------------------------------------------------;
        .code
        assume  cs:@code,ds:@data,es:nothing,ss:nothing
;-------------------------------------------------------;
;       enc     encode data                             ;
;-------------------------------------------------------;
r0      equ     al                      ;define remainder bytes
r1      equ     cl
r2      equ     dl
;
enc     proc    near
        push    bp
        mov     bp,sp
        push    ds
        push    es
        push    bx
        push    cx
        push    dx
        push    si
        push    di
        mov     ax,ss                   ;ds = ss
        mov     ds,ax
        les     di,[bp+4]
        mov     bx,mxc0
        mov     si,1024
;
enc0:   xor     ax,ax
        xor     cx,cx
        xor     dx,dx
encd    =       0
        xor     r1,es:[di+encd*1024]    ;do 1st 2 bytes
        mov     bl,r1
        xor     r0,[bx]
        xor     r2,[bx]
encd    =       encd+1
        xor     r0,es:[di+encd*1024]
        mov     bl,r0
        xor     r2,[bx]
        xor     r1,[bx]
encd    =       encd+1

        rept    9                       ;do remaining 27 bytes
        xor     r2,es:[di+encd*1024]
        mov     bl,r2
        xor     r1,[bx]
        xor     r0,[bx]
encd    =       encd+1
        xor     r1,es:[di+encd*1024]
        mov     bl,r1
        xor     r0,[bx]
        xor     r2,[bx]
encd    =       encd+1
        xor     r0,es:[di+encd*1024]
        mov     bl,r0
        xor     r2,[bx]
        xor     r1,[bx]
encd    =       encd+1
        endm

        mov     es:[di+29*1024],r2      ;store encoded data
        mov     es:[di+30*1024],r1
        mov     es:[di+31*1024],r0
        inc     di                      ;adv adr
        dec     si
        jnz     enc0
        pop     di
        pop     si
        pop     dx
        pop     cx
        pop     bx
        pop     es
        pop     ds
        pop     bp
        ret
enc     endp
;-------------------------------------------------------;
;       tbli    initialize table                        ;
;-------------------------------------------------------;
tbli    proc    near
        push    ds
        mov     ax,ss                   ;ds = ss
        mov     ds,ax
        mov     bx,exp2                 ;init exp2 table
        mov     al,01                   ; 00200h entries
        jmp     short tbli1             ; for faster mpy
tbli0:  add     al,al
        jnc     tbli1
        xor     al,87h
tbli1:  mov     [bx],al
        inc     bx
        cmp     bx,log2
        jne     tbli0
;
        mov     si,exp2                 ;init log2 table
        mov     byte ptr [bx],0ffh
        xor     ax,ax
tbli2:  mov     bl,[si]
        mov     [bx],al
        inc     si
        inc     ax
        cmp     ax,000ffh
        jne     tbli2
;
        mov     si,mxc0                 ;init mpy by c0 table
        mov     cx,0c000h
tbli3:  call    mpy
        mov     [si],al
        inc     si
        inc     cl
        jnz     tbli3
        pop     ds
        ret
tbli    endp
;-------------------------------------------------------;
;       mpy     al = ch*cl                              ;
;-------------------------------------------------------;
mpy     proc    near
        xor     ax,ax                   ;return 0 if 0
        or      cl,cl
        je      mpy0
        or      ch,ch
        je      mpy0
        mov     bx,log2                 ;bx=log[ch]+log[cl]
        mov     bl,ch
        mov     al,[bx]
        mov     bl,cl
        mov     bl,[bx]
        add     bx,ax
        mov     al,[bx+exp2-log2]       ;ax = ch*cl
mpy0:   ret
mpy     endp
;-------------------------------------------------------;
;       main                                            ;
;-------------------------------------------------------;
main    proc    far
        mov     ax,@data                ;set ds
        mov     ds,ax
        mov     ax,04a00h               ;free memory
        mov     bx,ss
        add     bx,stksz/16
        int     21h
        mov     ax,04800h               ;bx = # free paras
        mov     bx,0ffffh
        int     21h
        mov     ax,4800h                ;allocate
        int     21h
        cli
        mov     ss,ax                   ;ss = working segment
        mov     sp,wrksz
        sti
        add     ax,wrksz/16
        mov     es,ax                   ;es = bfr segment
        call    tbli                    ;init tables
        mov     di,bfr                  ;init bfr
        mov     cx,32768
        xor     ax,ax
        rep     stosb
        mov     di,bfr
        mov     cx,1024
        mov     al,01h
        rep     stosb
        xor     di,di
        push    es
        push    di
        call    enc                     ;encode
        add     sp,4
        mov     ax,4c00h
        int     21h
main    endp
        end     main

