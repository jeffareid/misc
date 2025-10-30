        page    80,128
        title   enc3
;-------------------------------------------------------;
;       enc3    386 enc code                            ;
;                                                       ;
;       Jeff Reid       11AUG89 17:15                   ;
;-------------------------------------------------------;
        .386p                   ;enable instructions
        .model flat,c           ;use C naming convention
;       include C libraries
        includelib      msvcrtd
;       includelib      oldnames
;       includelib      legacy_stdio_definitions.lib    ;for scanf, printf, ...

;-------------------------------------------------------;
;       unitialized data section                        ;
;-------------------------------------------------------;
        .data?
mxc0    word    65536 dup (?)           ;mpy by c0 table
bfr     byte    32768 dup (?)           ;bfr
exp2    byte    512 dup(?)              ;exp2 table
log2    byte    256 dup(?)              ;log2 table
;-------------------------------------------------------;
;       initialized data section                        ;
;-------------------------------------------------------;
        .data
;-------------------------------------------------------;
;       code section                                    ;
;-------------------------------------------------------;
        .code
;-------------------------------------------------------;
;       mpy     al = ch*cl                              ;
;-------------------------------------------------------;
mpy     proc    near
        xor     eax,eax                 ;return 0 if 0
        or      cl,cl
        je      mpy0
        or      ch,ch
        je      mpy0
        xor     ebx,ebx
        mov     bl,ch
        mov     al,log2[ebx]
        mov     bl,cl
        mov     bl,log2[ebx]
        add     ebx,eax
        mov     al,exp2[ebx]
mpy0:   ret
mpy     endp
;-------------------------------------------------------;
;       tbli    initialize table                        ;
;-------------------------------------------------------;
tbli    proc
        xor     ebx,ebx                 ;init exp2 table
        mov     al,001h                 ; 512 entries
        jmp     short tbli1             ; for faster mpy
tbli0:  add     al,al
        jnc     tbli1
        xor     al,87h
tbli1:  mov     exp2[ebx],al
        inc     ebx
        cmp     ebx,512
        jne     tbli0
;
        xor     ebx,ebx                 ;init log2 table
        xor     esi,esi
        mov     byte ptr log2[ebx],0ffh
        xor     eax,eax
tbli2:  mov     bl,exp2[esi]
        mov     log2[ebx],al
        inc     esi
        inc     eax
        cmp     eax,000ffh
        jne     tbli2
;
        lea     edi,mxc0                ;init mpy by c0 table
        mov     ecx,0c000h              ;ch = c0
        xor     edx,edx                 ;dh = row, dl = col
tbli3:  mov     cl,dl                   ;cl = col
        call    mpy
        stosb
        mov     cl,dh                   ;cl = row
        call    mpy
        stosb
        inc     dl
        jnz     tbli3                   ;loop till row done
        inc     dh
        jnz     tbli3                   ;loop till table done
        ret
tbli    endp
;-------------------------------------------------------;
;       enc     encode data                             ;
;                                                       ;
;       a single instruction multplies 2 bytes by c0    ;
;       and xor to working parity register such as      ;
;       xor     r2w,[ebx+r0d*2]                         ;
;-------------------------------------------------------;
r0w     equ     ax
r1w     equ     cx
r2w     equ     dx
r0d     equ     eax
r1d     equ     ecx
r2d     equ     edx
;
enc     proc    near
        lea     ebx,mxc0                ;mpy by c0
        mov     si,512                  ;# loops
;
enc0:   xor     eax,eax                 ;zero regs
        xor     ecx,ecx
        xor     edx,edx
encd    =       0
        xor     r1w,[edi+encd*1024]     ;do 1st 2 words
        xor     r0w,[ebx+r1d*2]
        xor     r2w,[ebx+r1d*2]
encd    =       encd+1
        xor     r0w,[edi+encd*1024]
        xor     r2w,[ebx+r0d*2]
        xor     r1w,[ebx+r0d*2]
encd    =       encd+1
        rept    9                       ;do remaining 27 words
        xor     r2w,[edi+encd*1024]
        xor     r1w,[ebx+r2d*2]
        xor     r0w,[ebx+r2d*2]
encd    =       encd+1
        xor     r1w,[edi+encd*1024]
        xor     r0w,[ebx+r1d*2]
        xor     r2w,[ebx+r1d*2]
encd    =       encd+1
        xor     r0w,[edi+encd*1024]
        xor     r2w,[ebx+r0d*2]
        xor     r1w,[ebx+r0d*2]
encd    =       encd+1
        endm
        mov     [edi+29*1024],r2w       ;store parity words
        mov     [edi+30*1024],r1w
        mov     [edi+31*1024],r0w
        add     edi,2                   ;adv adr
        dec     esi                     ;loop
        jnz     enc0
        ret
enc     endp
;-------------------------------------------------------;
;       main                                            ;
;-------------------------------------------------------;
main    proc
        push    ebp
        push    edi
        push    esi
        push    ebx
        call    tbli
        xor     eax,eax
        lea     edi,bfr
        mov     ecx,000002000h
        rep     stosd
        mov     eax,03040102h
        lea     edi,bfr
        mov     ecx,000000100h
        rep     stosd
        lea     edi,bfr
        call    enc
        lea     edi,bfr
        mov     eax,001010101h
        mov     ebx,29
main0:  mov     ecx,000000100h
        rep     stosd
        add     eax,001010101h
        dec     ebx
        jnz     main0
        lea     edi,bfr
        call    enc
        pop     ebx
        pop     esi
        pop     edi
        pop     ebp
        xor     eax,eax
        ret
main    endp

        end
 