
#include "basemul_core.inc"
#include "permute.inc"


.text

.if 0

a0 + w a1 + w^2 a2
a0 + w^2 a1 + w a2
=
a0 - a2 + (a1 - a2) w
a0 - a1 + (a2 - a1) w

.endif

.global __asm_3x2
.global ___asm_3x2
__asm_3x2:
___asm_3x2:

        mov  $17, %rax

        __asm_3x2_loop:

        vmovdqu 0x000(%rdi),  %ymm0
        vmovdqu 0x020(%rdi),  %ymm1
        vmovdqu 0x040(%rdi),  %ymm2

        vmovdqu 0x000(%rsi),  %ymm3
        vmovdqu 0x020(%rsi),  %ymm4
        vmovdqu 0x040(%rsi),  %ymm5

        vmovdqu 0x000(%rcx),  %ymm6
        vmovdqu 0x040(%rcx),  %ymm7
        vmovdqu 0x020(%rcx),  %ymm8
        vmovdqu 0x000(%rdx),  %ymm9

        barrett_reducex3 0, 1, 2, 0, 1, 2, 6, 7, 10, 11, 12
        barrett_reducex3 3, 4, 5, 3, 4, 5, 6, 7, 13, 14, 15

        vpaddw        %ymm3,  %ymm0, %ymm10
        vpsubw        %ymm3,  %ymm0, %ymm13
        vpaddw        %ymm4,  %ymm1, %ymm11
        vpsubw        %ymm4,  %ymm1, %ymm14
        vpaddw        %ymm5,  %ymm2, %ymm12
        vpsubw        %ymm5,  %ymm2, %ymm15

        vpaddw       %ymm11, %ymm10, %ymm0
        vpaddw       %ymm12,  %ymm0, %ymm0
        vpaddw       %ymm14, %ymm13, %ymm3
        vpaddw       %ymm15,  %ymm3, %ymm3

        barrett_reducex2 0, 3, 0, 3, 6, 7, 1, 4

        vmovdqu       %ymm0, 0x000(%rdi)
        vmovdqu       %ymm3, 0x000(%rsi)

        vpsubw       %ymm12, %ymm10,  %ymm1
        vpsubw       %ymm15, %ymm13,  %ymm4
        vpsubw       %ymm11, %ymm10,  %ymm2
        vpsubw       %ymm14, %ymm13,  %ymm5
        vpsubw       %ymm12, %ymm11, %ymm11
        vpsubw       %ymm15, %ymm14, %ymm14

        montgomery_mul 11, 11, 9, 6, 8, 0, 12
        montgomery_mul 14, 14, 9, 6, 8, 3, 15

        vpaddw       %ymm11,  %ymm1,  %ymm1
        vpaddw       %ymm14,  %ymm4,  %ymm4
        vpsubw       %ymm11,  %ymm2,  %ymm2
        vpsubw       %ymm14,  %ymm5,  %ymm5

        barrett_reducex2 1, 2, 1, 2, 6, 7, 14, 15
        barrett_reducex2 4, 5, 4, 5, 6, 7, 14, 15

        vmovdqu       %ymm1, 0x020(%rdi)
        vmovdqu       %ymm4, 0x020(%rsi)
        vmovdqu       %ymm2, 0x040(%rdi)
        vmovdqu       %ymm5, 0x040(%rsi)

        add          $0x060, %rdi
        add          $0x060, %rsi

        dec          %rax
        cmp          $0, %rax
        jne          __asm_3x2_loop


        ret

.global __asm_3x2_pre
.global ___asm_3x2_pre
__asm_3x2_pre:
___asm_3x2_pre:

        mov  $17, %rax

        __asm_3x2_pre_loop:

        vmovdqu 0x000(%rdi),  %ymm0
        vmovdqu 0x020(%rdi),  %ymm1
        vmovdqu 0x040(%rdi),  %ymm2

        vmovdqu 0x000(%rsi),  %ymm3
        vmovdqu 0x020(%rsi),  %ymm4
        vmovdqu 0x040(%rsi),  %ymm5

        vmovdqu 0x000(%rcx),  %ymm6
        vmovdqu 0x040(%rcx),  %ymm7
        vmovdqu 0x020(%rcx),  %ymm8
        vmovdqu 0x000(%rdx),  %ymm9

        barrett_reducex3 0, 1, 2, 0, 1, 2, 6, 7, 10, 11, 12
        barrett_reducex3 3, 4, 5, 3, 4, 5, 6, 7, 13, 14, 15

        vpaddw        %ymm3,  %ymm0, %ymm10
        vpsubw        %ymm3,  %ymm0, %ymm13
        vpaddw        %ymm4,  %ymm1, %ymm11
        vpsubw        %ymm4,  %ymm1, %ymm14
        vpaddw        %ymm5,  %ymm2, %ymm12
        vpsubw        %ymm5,  %ymm2, %ymm15

        vpaddw       %ymm11, %ymm10, %ymm0
        vpaddw       %ymm12,  %ymm0, %ymm0
        vpaddw       %ymm14, %ymm13, %ymm3
        vpaddw       %ymm15,  %ymm3, %ymm3

        barrett_reduce 3, 3, 6, 7, 4
        //barrett_reducex2 0, 3, 0, 3, 6, 7, 1, 4

        vmovdqu       %ymm0, 0x000(%rdi)
        vmovdqu       %ymm3, 0x000(%rsi)

        vpsubw       %ymm12, %ymm10,  %ymm1
        vpsubw       %ymm15, %ymm13,  %ymm4
        vpsubw       %ymm11, %ymm10,  %ymm2
        vpsubw       %ymm14, %ymm13,  %ymm5
        vpsubw       %ymm12, %ymm11, %ymm11
        vpsubw       %ymm15, %ymm14, %ymm14

        montgomery_mul 11, 11, 9, 6, 8, 0, 12
        montgomery_mul 14, 14, 9, 6, 8, 3, 15

        vpaddw       %ymm11,  %ymm1,  %ymm1
        vpaddw       %ymm14,  %ymm4,  %ymm4
        vpsubw       %ymm11,  %ymm2,  %ymm2
        vpsubw       %ymm14,  %ymm5,  %ymm5

        //barrett_reducex2 1, 2, 1, 2, 6, 7, 14, 15
        barrett_reducex2 4, 5, 4, 5, 6, 7, 14, 15

        vmovdqu       %ymm1, 0x020(%rdi)
        vmovdqu       %ymm4, 0x020(%rsi)
        vmovdqu       %ymm2, 0x040(%rdi)
        vmovdqu       %ymm5, 0x040(%rsi)

        add          $0x060, %rdi
        add          $0x060, %rsi

        dec          %rax
        cmp          $0, %rax
        jne          __asm_3x2_pre_loop


        ret

.global __asm_3x2_post
.global ___asm_3x2_post
__asm_3x2_post:
___asm_3x2_post:

        mov  $17, %rax

        __asm_3x2_post_loop:

        vmovdqu 0x000(%rdi),  %ymm0
        vmovdqu 0x020(%rdi),  %ymm1
        vmovdqu 0x040(%rdi),  %ymm2

        vmovdqu 0x000(%rsi),  %ymm3
        vmovdqu 0x020(%rsi),  %ymm4
        vmovdqu 0x040(%rsi),  %ymm5

        vmovdqu 0x000(%rcx),  %ymm6
        vmovdqu 0x040(%rcx),  %ymm7
        vmovdqu 0x020(%rcx),  %ymm8
        vmovdqu 0x000(%rdx),  %ymm9

        //barrett_reducex3 0, 1, 2, 0, 1, 2, 6, 7, 10, 11, 12
        //barrett_reducex3 3, 4, 5, 3, 4, 5, 6, 7, 13, 14, 15

        vpaddw        %ymm3,  %ymm0, %ymm10
        vpsubw        %ymm3,  %ymm0, %ymm13
        vpaddw        %ymm4,  %ymm1, %ymm11
        vpsubw        %ymm4,  %ymm1, %ymm14
        vpaddw        %ymm5,  %ymm2, %ymm12
        vpsubw        %ymm5,  %ymm2, %ymm15

        vpaddw       %ymm11, %ymm10, %ymm0
        vpaddw       %ymm12,  %ymm0, %ymm0
        vpaddw       %ymm14, %ymm13, %ymm3
        vpaddw       %ymm15,  %ymm3, %ymm3

        barrett_reducex2 0, 3, 0, 3, 6, 7, 1, 4

        vmovdqu       %ymm0, 0x000(%rdi)
        vmovdqu       %ymm3, 0x000(%rsi)

        vpsubw       %ymm12, %ymm10,  %ymm1
        vpsubw       %ymm15, %ymm13,  %ymm4
        vpsubw       %ymm11, %ymm10,  %ymm2
        vpsubw       %ymm14, %ymm13,  %ymm5
        vpsubw       %ymm12, %ymm11, %ymm11
        vpsubw       %ymm15, %ymm14, %ymm14

        montgomery_mul 11, 11, 9, 6, 8, 0, 12
        montgomery_mul 14, 14, 9, 6, 8, 3, 15

        vpaddw       %ymm11,  %ymm1,  %ymm1
        vpaddw       %ymm14,  %ymm4,  %ymm4
        vpsubw       %ymm11,  %ymm2,  %ymm2
        vpsubw       %ymm14,  %ymm5,  %ymm5

        barrett_reducex2 1, 2, 1, 2, 6, 7, 14, 15
        barrett_reducex2 4, 5, 4, 5, 6, 7, 14, 15

        vmovdqu       %ymm1, 0x020(%rdi)
        vmovdqu       %ymm4, 0x020(%rsi)
        vmovdqu       %ymm2, 0x040(%rdi)
        vmovdqu       %ymm5, 0x040(%rsi)

        add          $0x060, %rdi
        add          $0x060, %rsi

        dec          %rax
        cmp          $0, %rax
        jne          __asm_3x2_post_loop


        ret

.if 0

result of x^8 - 1

0x000
0x020
0x040
0x060
0x080
0x0a0
0x0c0
0x0e0

buffer for one of the polynomial in x^8 + 1

0x100
0x120
0x140
0x160
0x180
0x1a0
0x1c0
0x1e0

result of x^8 + 1

0x200
0x220
0x240
0x260
0x280
0x2a0
0x2c0
0x2e0

buffer for karatsuba negacyclic 8 x 8

0x300
0x320
0x340
0x360
0x380
0x3a0
0x3c0
0x3e0

0x400
0x420
0x440
0x460
0x480
0x4a0
0x4c0
0x4e0

0x500
0x520
0x540
0x560
0x580
0x5a0
0x5c0
0x5e0

.endif

.global __asm_rader17
.global ___asm_rader17
__asm_rader17:
___asm_rader17:

        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x700, %rsp

        vmovdqu     0(%rsi),  %ymm0
        vmovdqu       %ymm0, 0x600(%rsp)
        add     $0x020, %rdx
        add     $0x020, %rcx

.if 0

// g = 3
// map x to g^x
size_t rader_in_permute[17] = {
0,
3, 9, 10, 13,
5, 15, 11, 16,
14, 8, 7, 4,
12, 2, 6, 1
};

.endif

// 0, 96, 288, 320, 416, 160, 480, 352, 512, 448, 256, 224, 128, 384, 64, 192, 32

        vmovdqu    96(%rsi),  %ymm0
        vmovdqu   288(%rsi),  %ymm1
        vmovdqu   320(%rsi),  %ymm2
        vmovdqu   416(%rsi),  %ymm3
        vmovdqu   160(%rsi),  %ymm4
        vmovdqu   480(%rsi),  %ymm5
        vmovdqu   352(%rsi),  %ymm6
        vmovdqu   512(%rsi),  %ymm7

        vpsubw    448(%rsi),  %ymm0,  %ymm8
        vpaddw    448(%rsi),  %ymm0,  %ymm0
        vpsubw    256(%rsi),  %ymm1,  %ymm9
        vpaddw    256(%rsi),  %ymm1,  %ymm1
        vpsubw    224(%rsi),  %ymm2, %ymm10
        vpaddw    224(%rsi),  %ymm2,  %ymm2
        vpsubw    128(%rsi),  %ymm3, %ymm11
        vpaddw    128(%rsi),  %ymm3,  %ymm3
        vpsubw    384(%rsi),  %ymm4, %ymm12
        vpaddw    384(%rsi),  %ymm4,  %ymm4
        vpsubw     64(%rsi),  %ymm5, %ymm13
        vpaddw     64(%rsi),  %ymm5,  %ymm5
        vpsubw    192(%rsi),  %ymm6, %ymm14
        vpaddw    192(%rsi),  %ymm6,  %ymm6
        vpsubw     32(%rsi),  %ymm7, %ymm15
        vpaddw     32(%rsi),  %ymm7,  %ymm7

        vmovdqu       %ymm8, 0x100(%rsp)
        vmovdqu       %ymm9, 0x120(%rsp)
        vmovdqu      %ymm10, 0x140(%rsp)
        vmovdqu      %ymm11, 0x160(%rsp)
        vmovdqu      %ymm12, 0x180(%rsp)
        vmovdqu      %ymm13, 0x1a0(%rsp)
        vmovdqu      %ymm14, 0x1c0(%rsp)
        vmovdqu      %ymm15, 0x1e0(%rsp)

        vpaddw        %ymm4,  %ymm0,  %ymm8
        vpsubw        %ymm4,  %ymm0, %ymm12
        vpaddw        %ymm5,  %ymm1,  %ymm9
        vpsubw        %ymm5,  %ymm1, %ymm13
        vpaddw        %ymm6,  %ymm2, %ymm10
        vpsubw        %ymm6,  %ymm2, %ymm14
        vpaddw        %ymm7,  %ymm3, %ymm11
        vpsubw        %ymm7,  %ymm3, %ymm15

        vmovdqu      %ymm12, 0x080(%rsp)
        vmovdqu      %ymm13, 0x0a0(%rsp)
        vmovdqu      %ymm14, 0x0c0(%rsp)
        vmovdqu      %ymm15, 0x0e0(%rsp)

        vmovdqu 0x000( %r8), %ymm12
        vmovdqu 0x020( %r8), %ymm13

// ================================
// x^4 - 1

        vmovdqu  0x040( %r8),  %ymm13
        barrett_reducex2 0, 1,  8,  9, 12, 13, 14, 15
        barrett_reducex2 2, 3, 10, 11, 12, 13, 14, 15

        vpsubw        %ymm2,  %ymm0, %ymm10
        vpaddw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm1, %ymm11
        vpaddw        %ymm3,  %ymm1,  %ymm1
        vpaddw        %ymm1,  %ymm0,  %ymm8
        vpsubw        %ymm1,  %ymm0,  %ymm9

        vmovdqu       %ymm8, 0x620(%rsp)

        vmovdqu 0x000(%rdx),  %ymm4
        vmovdqu 0x020(%rdx),  %ymm5
        vmovdqu 0x040(%rdx),  %ymm6
        vmovdqu 0x060(%rdx),  %ymm7

        vpmullw 0x000(%rcx),  %ymm8, %ymm14
        vpmullw 0x020(%rcx),  %ymm9, %ymm15
        vpmulhw       %ymm4,  %ymm8,  %ymm4
        vpmulhw       %ymm5,  %ymm9,  %ymm5
        vpmulhw      %ymm14, %ymm12, %ymm14
        vpmulhw      %ymm15, %ymm12, %ymm15
        vpsubw       %ymm14,  %ymm4,  %ymm4
        vpsubw       %ymm15,  %ymm5,  %ymm5

        vpmullw 0x040(%rcx), %ymm10, %ymm14
        vpmullw 0x060(%rcx), %ymm10, %ymm15
        vpmullw 0x040(%rcx), %ymm11,  %ymm8
        vpmullw 0x060(%rcx), %ymm11,  %ymm9
        vpmulhw       %ymm6, %ymm10,  %ymm0
        vpmulhw       %ymm7, %ymm10,  %ymm1
        vpmulhw       %ymm6, %ymm11,  %ymm2
        vpmulhw       %ymm7, %ymm11,  %ymm3
        vpmulhw      %ymm14, %ymm12, %ymm14
        vpmulhw      %ymm15, %ymm12, %ymm15
        vpmulhw       %ymm8, %ymm12,  %ymm8
        vpmulhw       %ymm9, %ymm12,  %ymm9
        vpsubw       %ymm14,  %ymm0,  %ymm0
        vpsubw       %ymm15,  %ymm1,  %ymm1
        vpsubw        %ymm8,  %ymm2,  %ymm2
        vpsubw        %ymm9,  %ymm3,  %ymm3

        vpsubw        %ymm3,  %ymm0, %ymm6
        vpaddw        %ymm2,  %ymm1, %ymm7

        vpaddw        %ymm5,  %ymm4,  %ymm0
        vpsubw        %ymm5,  %ymm4,  %ymm1
        vpsubw        %ymm6,  %ymm0,  %ymm2
        vpaddw        %ymm6,  %ymm0,  %ymm0
        vpsubw        %ymm7,  %ymm1,  %ymm3
        vpaddw        %ymm7,  %ymm1,  %ymm1

        vmovdqu  0x040( %r8),  %ymm13

        barrett_reducex2 0, 1, 0, 1, 12, 13, 14, 15
        barrett_reducex2 2, 3, 2, 3, 12, 13, 14, 15

        vmovdqu       %ymm0, 0x000(%rsp)
        vmovdqu       %ymm1, 0x020(%rsp)
        vmovdqu       %ymm2, 0x040(%rsp)
        vmovdqu       %ymm3, 0x060(%rsp)

        vmovdqu  0x020( %r8),  %ymm13

// ================================
// x^4 + 1

        vmovdqu 0x080(%rsp),  %ymm8
        vmovdqu 0x0a0(%rsp),  %ymm9
        vmovdqu 0x0c0(%rsp), %ymm10
        vmovdqu 0x0e0(%rsp), %ymm11

        vpmullw 0x080(%rcx),  %ymm8,  %ymm0
        vpmullw 0x0a0(%rcx), %ymm11,  %ymm1
        vpmullw 0x0c0(%rcx), %ymm10,  %ymm2
        vpmullw 0x0e0(%rcx),  %ymm9,  %ymm3
        vpmulhw       %ymm0, %ymm12,  %ymm0
        vpmulhw       %ymm1, %ymm12,  %ymm1
        vpmulhw       %ymm2, %ymm12,  %ymm2
        vpmulhw       %ymm3, %ymm12,  %ymm3
        vpmulhw 0x080(%rdx),  %ymm8, %ymm12
        vpmulhw 0x0a0(%rdx), %ymm11, %ymm13
        vpmulhw 0x0c0(%rdx), %ymm10, %ymm14
        vpmulhw 0x0e0(%rdx),  %ymm9, %ymm15
        vpsubw        %ymm0, %ymm12,  %ymm0
        vpsubw        %ymm1, %ymm13,  %ymm1
        vpsubw        %ymm2, %ymm14,  %ymm2
        vpsubw        %ymm3, %ymm15,  %ymm3

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpsubw        %ymm1,  %ymm0,  %ymm0
        vpsubw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm0,  %ymm0

        barrett_reduce 0, 0, 12, 13, 14
        vmovdqu       %ymm0, 0x080(%rsp)

        vpmullw 0x080(%rcx),  %ymm9,  %ymm4
        vpmullw 0x0a0(%rcx),  %ymm8,  %ymm5
        vpmullw 0x0c0(%rcx), %ymm11,  %ymm6
        vpmullw 0x0e0(%rcx), %ymm10,  %ymm7
        vpmulhw       %ymm4, %ymm12,  %ymm4
        vpmulhw       %ymm5, %ymm12,  %ymm5
        vpmulhw       %ymm6, %ymm12,  %ymm6
        vpmulhw       %ymm7, %ymm12,  %ymm7
        vpmulhw 0x080(%rdx),  %ymm9, %ymm12
        vpmulhw 0x0a0(%rdx),  %ymm8, %ymm13
        vpmulhw 0x0c0(%rdx), %ymm11, %ymm14
        vpmulhw 0x0e0(%rdx), %ymm10, %ymm15
        vpsubw        %ymm4, %ymm12,  %ymm4
        vpsubw        %ymm5, %ymm13,  %ymm5
        vpsubw        %ymm6, %ymm14,  %ymm6
        vpsubw        %ymm7, %ymm15,  %ymm7

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm5,  %ymm4,  %ymm4
        vpsubw        %ymm6,  %ymm4,  %ymm4
        vpsubw        %ymm7,  %ymm4,  %ymm4

        barrett_reduce 4, 4, 12, 13, 14
        vmovdqu       %ymm4, 0x0a0(%rsp)

        vpmullw 0x080(%rcx), %ymm10,  %ymm0
        vpmullw 0x0a0(%rcx),  %ymm9,  %ymm1
        vpmullw 0x0c0(%rcx),  %ymm8,  %ymm2
        vpmullw 0x0e0(%rcx), %ymm11,  %ymm3
        vpmulhw       %ymm0, %ymm12,  %ymm0
        vpmulhw       %ymm1, %ymm12,  %ymm1
        vpmulhw       %ymm2, %ymm12,  %ymm2
        vpmulhw       %ymm3, %ymm12,  %ymm3
        vpmulhw 0x080(%rdx), %ymm10, %ymm12
        vpmulhw 0x0a0(%rdx),  %ymm9, %ymm13
        vpmulhw 0x0c0(%rdx),  %ymm8, %ymm14
        vpmulhw 0x0e0(%rdx), %ymm11, %ymm15
        vpsubw        %ymm0, %ymm12,  %ymm0
        vpsubw        %ymm1, %ymm13,  %ymm1
        vpsubw        %ymm2, %ymm14,  %ymm2
        vpsubw        %ymm3, %ymm15,  %ymm3

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm1,  %ymm0,  %ymm0
        vpaddw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm0,  %ymm0

        barrett_reduce 0, 0, 12, 13, 14
        vmovdqu       %ymm0, 0x0c0(%rsp)

        vpmullw 0x080(%rcx), %ymm11,  %ymm4
        vpmullw 0x0a0(%rcx), %ymm10,  %ymm5
        vpmullw 0x0c0(%rcx),  %ymm9,  %ymm6
        vpmullw 0x0e0(%rcx),  %ymm8,  %ymm7
        vpmulhw       %ymm4, %ymm12,  %ymm4
        vpmulhw       %ymm5, %ymm12,  %ymm5
        vpmulhw       %ymm6, %ymm12,  %ymm6
        vpmulhw       %ymm7, %ymm12,  %ymm7
        vpmulhw 0x080(%rdx), %ymm11, %ymm12
        vpmulhw 0x0a0(%rdx), %ymm10, %ymm13
        vpmulhw 0x0c0(%rdx),  %ymm9, %ymm14
        vpmulhw 0x0e0(%rdx),  %ymm8, %ymm15
        vpsubw        %ymm4, %ymm12,  %ymm4
        vpsubw        %ymm5, %ymm13,  %ymm5
        vpsubw        %ymm6, %ymm14,  %ymm6
        vpsubw        %ymm7, %ymm15,  %ymm7

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm5,  %ymm4,  %ymm4
        vpaddw        %ymm6,  %ymm4,  %ymm4
        vpaddw        %ymm7,  %ymm4,  %ymm4

        barrett_reduce 4, 4, 12, 13, 14
        vmovdqu       %ymm4, 0x0e0(%rsp)

// ================================
// x^8 + 1

        __asm_negacyclic_karatsuba8_precompute_core \
                %rsp, %rsp, %rdx, %rcx, %rsp, 0x200, 0x100, 0x100, 0x100, 0x300, \
                0x000( %r8), 0x020( %r8), 0x040( %r8)

        vmovdqu 0x000(%rsp),  %ymm0
        vmovdqu 0x020(%rsp),  %ymm1
        vmovdqu 0x040(%rsp),  %ymm2
        vmovdqu 0x060(%rsp),  %ymm3

        vpsubw  0x080(%rsp),  %ymm0,  %ymm4
        vpaddw  0x080(%rsp),  %ymm0,  %ymm0
        vpsubw  0x0a0(%rsp),  %ymm1,  %ymm5
        vpaddw  0x0a0(%rsp),  %ymm1,  %ymm1
        vpsubw  0x0c0(%rsp),  %ymm2,  %ymm6
        vpaddw  0x0c0(%rsp),  %ymm2,  %ymm2
        vpsubw  0x0e0(%rsp),  %ymm3,  %ymm7
        vpaddw  0x0e0(%rsp),  %ymm3,  %ymm3

        vpsubw  0x200(%rsp),  %ymm0,  %ymm8
        vpaddw  0x200(%rsp),  %ymm0,  %ymm0
        vpsubw  0x220(%rsp),  %ymm1,  %ymm9
        vpaddw  0x220(%rsp),  %ymm1,  %ymm1
        vpsubw  0x240(%rsp),  %ymm2, %ymm10
        vpaddw  0x240(%rsp),  %ymm2,  %ymm2
        vpsubw  0x260(%rsp),  %ymm3, %ymm11
        vpaddw  0x260(%rsp),  %ymm3,  %ymm3
        vpsubw  0x280(%rsp),  %ymm4, %ymm12
        vpaddw  0x280(%rsp),  %ymm4,  %ymm4
        vpsubw  0x2a0(%rsp),  %ymm5, %ymm13
        vpaddw  0x2a0(%rsp),  %ymm5,  %ymm5
        vpsubw  0x2c0(%rsp),  %ymm6, %ymm14
        vpaddw  0x2c0(%rsp),  %ymm6,  %ymm6
        vpsubw  0x2e0(%rsp),  %ymm7, %ymm15
        vpaddw  0x2e0(%rsp),  %ymm7,  %ymm7

        vpaddw  0x600(%rsp),  %ymm0,  %ymm0
        vpaddw  0x600(%rsp),  %ymm1,  %ymm1
        vpaddw  0x600(%rsp),  %ymm2,  %ymm2
        vpaddw  0x600(%rsp),  %ymm3,  %ymm3
        vpaddw  0x600(%rsp),  %ymm4,  %ymm4
        vpaddw  0x600(%rsp),  %ymm5,  %ymm5
        vpaddw  0x600(%rsp),  %ymm6,  %ymm6
        vpaddw  0x600(%rsp),  %ymm7,  %ymm7
        vpaddw  0x600(%rsp),  %ymm8,  %ymm8
        vpaddw  0x600(%rsp),  %ymm9,  %ymm9
        vpaddw  0x600(%rsp), %ymm10, %ymm10
        vpaddw  0x600(%rsp), %ymm11, %ymm11
        vpaddw  0x600(%rsp), %ymm12, %ymm12
        vpaddw  0x600(%rsp), %ymm13, %ymm13
        vpaddw  0x600(%rsp), %ymm14, %ymm14
        vpaddw  0x600(%rsp), %ymm15, %ymm15


.if 0

// g = 3
// map x to g^{16 - x}
size_t rader_out_permute[17] = {
0,
6, 2, 12, 4,
7, 8, 14, 16,
11, 15, 5, 13,
10, 9, 3, 1
};

000

1:
020
040
060
080

5:
0a0
0c0
0e0
100

9:
120
140
160
180

13:
1a0
1c0
1e0
200

.endif

// 0, 192, 64, 384, 128, 224, 256, 448, 512, 352, 480, 160, 416, 320, 288, 96, 32

        vmovdqu       %ymm0,   192(%rdi)
        vmovdqu       %ymm1,    64(%rdi)
        vmovdqu       %ymm2,   384(%rdi)
        vmovdqu       %ymm3,   128(%rdi)
        vmovdqu       %ymm4,   224(%rdi)
        vmovdqu       %ymm5,   256(%rdi)
        vmovdqu       %ymm6,   448(%rdi)
        vmovdqu       %ymm7,   512(%rdi)
        vmovdqu       %ymm8,   352(%rdi)
        vmovdqu       %ymm9,   480(%rdi)
        vmovdqu      %ymm10,   160(%rdi)
        vmovdqu      %ymm11,   416(%rdi)
        vmovdqu      %ymm12,   320(%rdi)
        vmovdqu      %ymm13,   288(%rdi)
        vmovdqu      %ymm14,    96(%rdi)
        vmovdqu      %ymm15,    32(%rdi)

        vmovdqu 0x620(%rsp),  %ymm0
        vpaddw  0x600(%rsp),  %ymm0,  %ymm0
        vmovdqu       %ymm0,      0(%rdi)



        vzeroupper
        leave
        ret

.global __asm_rader17_scalei
.global ___asm_rader17_scalei
__asm_rader17_scalei:
___asm_rader17_scalei:

        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x700, %rsp

        vmovdqu     0(%rsi),  %ymm0
        vmovdqu       %ymm0, 0x600(%rsp)
        add     $0x020, %rdx
        add     $0x020, %rcx

.if 0

// g = 3
// map x to g^x
size_t rader_in_permute[17] = {
0,
3, 9, 10, 13,
5, 15, 11, 16,
14, 8, 7, 4,
12, 2, 6, 1
};

.endif

// 0, 288, 864, 960, 1248, 480, 1440, 1056, 1536, 1344, 768, 672, 384, 1152, 192, 576, 96

        vmovdqu   288(%rsi),  %ymm0
        vmovdqu   864(%rsi),  %ymm1
        vmovdqu   960(%rsi),  %ymm2
        vmovdqu  1248(%rsi),  %ymm3
        vmovdqu   480(%rsi),  %ymm4
        vmovdqu  1440(%rsi),  %ymm5
        vmovdqu  1056(%rsi),  %ymm6
        vmovdqu  1536(%rsi),  %ymm7

        vpsubw   1344(%rsi),  %ymm0,  %ymm8
        vpaddw   1344(%rsi),  %ymm0,  %ymm0
        vpsubw    768(%rsi),  %ymm1,  %ymm9
        vpaddw    768(%rsi),  %ymm1,  %ymm1
        vpsubw    672(%rsi),  %ymm2, %ymm10
        vpaddw    672(%rsi),  %ymm2,  %ymm2
        vpsubw    384(%rsi),  %ymm3, %ymm11
        vpaddw    384(%rsi),  %ymm3,  %ymm3
        vpsubw   1152(%rsi),  %ymm4, %ymm12
        vpaddw   1152(%rsi),  %ymm4,  %ymm4
        vpsubw    192(%rsi),  %ymm5, %ymm13
        vpaddw    192(%rsi),  %ymm5,  %ymm5
        vpsubw    576(%rsi),  %ymm6, %ymm14
        vpaddw    576(%rsi),  %ymm6,  %ymm6
        vpsubw     96(%rsi),  %ymm7, %ymm15
        vpaddw     96(%rsi),  %ymm7,  %ymm7

        vmovdqu       %ymm8, 0x100(%rsp)
        vmovdqu       %ymm9, 0x120(%rsp)
        vmovdqu      %ymm10, 0x140(%rsp)
        vmovdqu      %ymm11, 0x160(%rsp)
        vmovdqu      %ymm12, 0x180(%rsp)
        vmovdqu      %ymm13, 0x1a0(%rsp)
        vmovdqu      %ymm14, 0x1c0(%rsp)
        vmovdqu      %ymm15, 0x1e0(%rsp)

        vpaddw        %ymm4,  %ymm0,  %ymm8
        vpsubw        %ymm4,  %ymm0, %ymm12
        vpaddw        %ymm5,  %ymm1,  %ymm9
        vpsubw        %ymm5,  %ymm1, %ymm13
        vpaddw        %ymm6,  %ymm2, %ymm10
        vpsubw        %ymm6,  %ymm2, %ymm14
        vpaddw        %ymm7,  %ymm3, %ymm11
        vpsubw        %ymm7,  %ymm3, %ymm15

        vmovdqu      %ymm12, 0x080(%rsp)
        vmovdqu      %ymm13, 0x0a0(%rsp)
        vmovdqu      %ymm14, 0x0c0(%rsp)
        vmovdqu      %ymm15, 0x0e0(%rsp)

        vmovdqu 0x000( %r8), %ymm12
        vmovdqu 0x020( %r8), %ymm13

// ================================
// x^4 - 1

        vmovdqu  0x040( %r8),  %ymm13
        barrett_reducex2 0, 1,  8,  9, 12, 13, 14, 15
        barrett_reducex2 2, 3, 10, 11, 12, 13, 14, 15

        vpsubw        %ymm2,  %ymm0, %ymm10
        vpaddw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm1, %ymm11
        vpaddw        %ymm3,  %ymm1,  %ymm1
        vpaddw        %ymm1,  %ymm0,  %ymm8
        vpsubw        %ymm1,  %ymm0,  %ymm9

        vmovdqu       %ymm8, 0x620(%rsp)

        vmovdqu 0x000(%rdx),  %ymm4
        vmovdqu 0x020(%rdx),  %ymm5
        vmovdqu 0x040(%rdx),  %ymm6
        vmovdqu 0x060(%rdx),  %ymm7

        vpmullw 0x000(%rcx),  %ymm8, %ymm14
        vpmullw 0x020(%rcx),  %ymm9, %ymm15
        vpmulhw       %ymm4,  %ymm8,  %ymm4
        vpmulhw       %ymm5,  %ymm9,  %ymm5
        vpmulhw      %ymm14, %ymm12, %ymm14
        vpmulhw      %ymm15, %ymm12, %ymm15
        vpsubw       %ymm14,  %ymm4,  %ymm4
        vpsubw       %ymm15,  %ymm5,  %ymm5

        vpmullw 0x040(%rcx), %ymm10, %ymm14
        vpmullw 0x060(%rcx), %ymm10, %ymm15
        vpmullw 0x040(%rcx), %ymm11,  %ymm8
        vpmullw 0x060(%rcx), %ymm11,  %ymm9
        vpmulhw       %ymm6, %ymm10,  %ymm0
        vpmulhw       %ymm7, %ymm10,  %ymm1
        vpmulhw       %ymm6, %ymm11,  %ymm2
        vpmulhw       %ymm7, %ymm11,  %ymm3
        vpmulhw      %ymm14, %ymm12, %ymm14
        vpmulhw      %ymm15, %ymm12, %ymm15
        vpmulhw       %ymm8, %ymm12,  %ymm8
        vpmulhw       %ymm9, %ymm12,  %ymm9
        vpsubw       %ymm14,  %ymm0,  %ymm0
        vpsubw       %ymm15,  %ymm1,  %ymm1
        vpsubw        %ymm8,  %ymm2,  %ymm2
        vpsubw        %ymm9,  %ymm3,  %ymm3

        vpsubw        %ymm3,  %ymm0, %ymm6
        vpaddw        %ymm2,  %ymm1, %ymm7

        vpaddw        %ymm5,  %ymm4,  %ymm0
        vpsubw        %ymm5,  %ymm4,  %ymm1
        vpsubw        %ymm6,  %ymm0,  %ymm2
        vpaddw        %ymm6,  %ymm0,  %ymm0
        vpsubw        %ymm7,  %ymm1,  %ymm3
        vpaddw        %ymm7,  %ymm1,  %ymm1

        vmovdqu  0x040( %r8),  %ymm13

        barrett_reducex2 0, 1, 0, 1, 12, 13, 14, 15
        barrett_reducex2 2, 3, 2, 3, 12, 13, 14, 15

        vmovdqu       %ymm0, 0x000(%rsp)
        vmovdqu       %ymm1, 0x020(%rsp)
        vmovdqu       %ymm2, 0x040(%rsp)
        vmovdqu       %ymm3, 0x060(%rsp)

        vmovdqu  0x020( %r8),  %ymm13

// ================================
// x^4 + 1

        vmovdqu 0x080(%rsp),  %ymm8
        vmovdqu 0x0a0(%rsp),  %ymm9
        vmovdqu 0x0c0(%rsp), %ymm10
        vmovdqu 0x0e0(%rsp), %ymm11

        vpmullw 0x080(%rcx),  %ymm8,  %ymm0
        vpmullw 0x0a0(%rcx), %ymm11,  %ymm1
        vpmullw 0x0c0(%rcx), %ymm10,  %ymm2
        vpmullw 0x0e0(%rcx),  %ymm9,  %ymm3
        vpmulhw       %ymm0, %ymm12,  %ymm0
        vpmulhw       %ymm1, %ymm12,  %ymm1
        vpmulhw       %ymm2, %ymm12,  %ymm2
        vpmulhw       %ymm3, %ymm12,  %ymm3
        vpmulhw 0x080(%rdx),  %ymm8, %ymm12
        vpmulhw 0x0a0(%rdx), %ymm11, %ymm13
        vpmulhw 0x0c0(%rdx), %ymm10, %ymm14
        vpmulhw 0x0e0(%rdx),  %ymm9, %ymm15
        vpsubw        %ymm0, %ymm12,  %ymm0
        vpsubw        %ymm1, %ymm13,  %ymm1
        vpsubw        %ymm2, %ymm14,  %ymm2
        vpsubw        %ymm3, %ymm15,  %ymm3

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpsubw        %ymm1,  %ymm0,  %ymm0
        vpsubw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm0,  %ymm0

        barrett_reduce 0, 0, 12, 13, 14
        vmovdqu       %ymm0, 0x080(%rsp)

        vpmullw 0x080(%rcx),  %ymm9,  %ymm4
        vpmullw 0x0a0(%rcx),  %ymm8,  %ymm5
        vpmullw 0x0c0(%rcx), %ymm11,  %ymm6
        vpmullw 0x0e0(%rcx), %ymm10,  %ymm7
        vpmulhw       %ymm4, %ymm12,  %ymm4
        vpmulhw       %ymm5, %ymm12,  %ymm5
        vpmulhw       %ymm6, %ymm12,  %ymm6
        vpmulhw       %ymm7, %ymm12,  %ymm7
        vpmulhw 0x080(%rdx),  %ymm9, %ymm12
        vpmulhw 0x0a0(%rdx),  %ymm8, %ymm13
        vpmulhw 0x0c0(%rdx), %ymm11, %ymm14
        vpmulhw 0x0e0(%rdx), %ymm10, %ymm15
        vpsubw        %ymm4, %ymm12,  %ymm4
        vpsubw        %ymm5, %ymm13,  %ymm5
        vpsubw        %ymm6, %ymm14,  %ymm6
        vpsubw        %ymm7, %ymm15,  %ymm7

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm5,  %ymm4,  %ymm4
        vpsubw        %ymm6,  %ymm4,  %ymm4
        vpsubw        %ymm7,  %ymm4,  %ymm4

        barrett_reduce 4, 4, 12, 13, 14
        vmovdqu       %ymm4, 0x0a0(%rsp)

        vpmullw 0x080(%rcx), %ymm10,  %ymm0
        vpmullw 0x0a0(%rcx),  %ymm9,  %ymm1
        vpmullw 0x0c0(%rcx),  %ymm8,  %ymm2
        vpmullw 0x0e0(%rcx), %ymm11,  %ymm3
        vpmulhw       %ymm0, %ymm12,  %ymm0
        vpmulhw       %ymm1, %ymm12,  %ymm1
        vpmulhw       %ymm2, %ymm12,  %ymm2
        vpmulhw       %ymm3, %ymm12,  %ymm3
        vpmulhw 0x080(%rdx), %ymm10, %ymm12
        vpmulhw 0x0a0(%rdx),  %ymm9, %ymm13
        vpmulhw 0x0c0(%rdx),  %ymm8, %ymm14
        vpmulhw 0x0e0(%rdx), %ymm11, %ymm15
        vpsubw        %ymm0, %ymm12,  %ymm0
        vpsubw        %ymm1, %ymm13,  %ymm1
        vpsubw        %ymm2, %ymm14,  %ymm2
        vpsubw        %ymm3, %ymm15,  %ymm3

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm1,  %ymm0,  %ymm0
        vpaddw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm0,  %ymm0

        barrett_reduce 0, 0, 12, 13, 14
        vmovdqu       %ymm0, 0x0c0(%rsp)

        vpmullw 0x080(%rcx), %ymm11,  %ymm4
        vpmullw 0x0a0(%rcx), %ymm10,  %ymm5
        vpmullw 0x0c0(%rcx),  %ymm9,  %ymm6
        vpmullw 0x0e0(%rcx),  %ymm8,  %ymm7
        vpmulhw       %ymm4, %ymm12,  %ymm4
        vpmulhw       %ymm5, %ymm12,  %ymm5
        vpmulhw       %ymm6, %ymm12,  %ymm6
        vpmulhw       %ymm7, %ymm12,  %ymm7
        vpmulhw 0x080(%rdx), %ymm11, %ymm12
        vpmulhw 0x0a0(%rdx), %ymm10, %ymm13
        vpmulhw 0x0c0(%rdx),  %ymm9, %ymm14
        vpmulhw 0x0e0(%rdx),  %ymm8, %ymm15
        vpsubw        %ymm4, %ymm12,  %ymm4
        vpsubw        %ymm5, %ymm13,  %ymm5
        vpsubw        %ymm6, %ymm14,  %ymm6
        vpsubw        %ymm7, %ymm15,  %ymm7

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm5,  %ymm4,  %ymm4
        vpaddw        %ymm6,  %ymm4,  %ymm4
        vpaddw        %ymm7,  %ymm4,  %ymm4

        barrett_reduce 4, 4, 12, 13, 14
        vmovdqu       %ymm4, 0x0e0(%rsp)

// ================================
// x^8 + 1

        __asm_negacyclic_karatsuba8_precompute_core \
                %rsp, %rsp, %rdx, %rcx, %rsp, 0x200, 0x100, 0x100, 0x100, 0x300, \
                0x000( %r8), 0x020( %r8), 0x040( %r8)

        vmovdqu 0x000(%rsp),  %ymm0
        vmovdqu 0x020(%rsp),  %ymm1
        vmovdqu 0x040(%rsp),  %ymm2
        vmovdqu 0x060(%rsp),  %ymm3

        vpsubw  0x080(%rsp),  %ymm0,  %ymm4
        vpaddw  0x080(%rsp),  %ymm0,  %ymm0
        vpsubw  0x0a0(%rsp),  %ymm1,  %ymm5
        vpaddw  0x0a0(%rsp),  %ymm1,  %ymm1
        vpsubw  0x0c0(%rsp),  %ymm2,  %ymm6
        vpaddw  0x0c0(%rsp),  %ymm2,  %ymm2
        vpsubw  0x0e0(%rsp),  %ymm3,  %ymm7
        vpaddw  0x0e0(%rsp),  %ymm3,  %ymm3

        vpsubw  0x200(%rsp),  %ymm0,  %ymm8
        vpaddw  0x200(%rsp),  %ymm0,  %ymm0
        vpsubw  0x220(%rsp),  %ymm1,  %ymm9
        vpaddw  0x220(%rsp),  %ymm1,  %ymm1
        vpsubw  0x240(%rsp),  %ymm2, %ymm10
        vpaddw  0x240(%rsp),  %ymm2,  %ymm2
        vpsubw  0x260(%rsp),  %ymm3, %ymm11
        vpaddw  0x260(%rsp),  %ymm3,  %ymm3
        vpsubw  0x280(%rsp),  %ymm4, %ymm12
        vpaddw  0x280(%rsp),  %ymm4,  %ymm4
        vpsubw  0x2a0(%rsp),  %ymm5, %ymm13
        vpaddw  0x2a0(%rsp),  %ymm5,  %ymm5
        vpsubw  0x2c0(%rsp),  %ymm6, %ymm14
        vpaddw  0x2c0(%rsp),  %ymm6,  %ymm6
        vpsubw  0x2e0(%rsp),  %ymm7, %ymm15
        vpaddw  0x2e0(%rsp),  %ymm7,  %ymm7

        vpaddw  0x600(%rsp),  %ymm0,  %ymm0
        vpaddw  0x600(%rsp),  %ymm1,  %ymm1
        vpaddw  0x600(%rsp),  %ymm2,  %ymm2
        vpaddw  0x600(%rsp),  %ymm3,  %ymm3
        vpaddw  0x600(%rsp),  %ymm4,  %ymm4
        vpaddw  0x600(%rsp),  %ymm5,  %ymm5
        vpaddw  0x600(%rsp),  %ymm6,  %ymm6
        vpaddw  0x600(%rsp),  %ymm7,  %ymm7
        vpaddw  0x600(%rsp),  %ymm8,  %ymm8
        vpaddw  0x600(%rsp),  %ymm9,  %ymm9
        vpaddw  0x600(%rsp), %ymm10, %ymm10
        vpaddw  0x600(%rsp), %ymm11, %ymm11
        vpaddw  0x600(%rsp), %ymm12, %ymm12
        vpaddw  0x600(%rsp), %ymm13, %ymm13
        vpaddw  0x600(%rsp), %ymm14, %ymm14
        vpaddw  0x600(%rsp), %ymm15, %ymm15


.if 0

// g = 3
// map x to g^{16 - x}
size_t rader_out_permute[17] = {
0,
6, 2, 12, 4,
7, 8, 14, 16,
11, 15, 5, 13,
10, 9, 3, 1
};

.endif

// 0, 576, 192, 1152, 384, 672, 768, 1344, 1536, 1056, 1440, 480, 1248, 960, 864, 288, 96

        vmovdqu       %ymm0,   576(%rdi)
        vmovdqu       %ymm1,   192(%rdi)
        vmovdqu       %ymm2,  1152(%rdi)
        vmovdqu       %ymm3,   384(%rdi)
        vmovdqu       %ymm4,   672(%rdi)
        vmovdqu       %ymm5,   768(%rdi)
        vmovdqu       %ymm6,  1344(%rdi)
        vmovdqu       %ymm7,  1536(%rdi)
        vmovdqu       %ymm8,  1056(%rdi)
        vmovdqu       %ymm9,  1440(%rdi)
        vmovdqu      %ymm10,   480(%rdi)
        vmovdqu      %ymm11,  1248(%rdi)
        vmovdqu      %ymm12,   960(%rdi)
        vmovdqu      %ymm13,   864(%rdi)
        vmovdqu      %ymm14,   288(%rdi)
        vmovdqu      %ymm15,    96(%rdi)

        vmovdqu 0x620(%rsp),  %ymm0
        vpaddw  0x600(%rsp),  %ymm0,  %ymm0
        vmovdqu       %ymm0,      0(%rdi)

        vzeroupper
        leave
        ret

.global __asm_rader17_pre_tbl
.global ___asm_rader17_pre_tbl
__asm_rader17_pre_tbl:
___asm_rader17_pre_tbl:

        push   %rbp
        push   %r12
        push   %r13
        push   %r14
        push   %r15
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x700, %rsp

        mov        0*8( %r9), %r12

        vmovdqu (%rsi, %r12),  %ymm0
        vmovdqu        %ymm0, 0x600(%rsp)
        add     $0x020, %rdx
        add     $0x020, %rcx

.if 0

// g = 3
// map x to g^x
size_t rader_in_permute[17] = {
0,
3, 9, 10, 13,
5, 15, 11, 16,
14, 8, 7, 4,
12, 2, 6, 1
};

.endif

// 0, 288, 864, 960, 1248, 480, 1440, 1056, 1536, 1344, 768, 672, 384, 1152, 192, 576, 96

        mov        3*8( %r9), %r12
        mov        9*8( %r9), %r13
        mov       10*8( %r9), %r14
        mov       13*8( %r9), %r15

        vmovdqu (%rsi, %r12),  %ymm0
        vmovdqu (%rsi, %r13),  %ymm1
        vmovdqu (%rsi, %r14),  %ymm2
        vmovdqu (%rsi, %r15),  %ymm3

        mov        5*8( %r9), %r12
        mov       15*8( %r9), %r13
        mov       11*8( %r9), %r14
        mov       16*8( %r9), %r15

        vmovdqu (%rsi, %r12),  %ymm4
        vmovdqu (%rsi, %r13),  %ymm5
        vmovdqu (%rsi, %r14),  %ymm6
        vmovdqu (%rsi, %r15),  %ymm7

        mov       14*8( %r9), %r12
        mov        8*8( %r9), %r13
        mov        7*8( %r9), %r14
        mov        4*8( %r9), %r15

        vpsubw (%rsi, %r12),  %ymm0,  %ymm8
        vpaddw (%rsi, %r12),  %ymm0,  %ymm0
        vpsubw (%rsi, %r13),  %ymm1,  %ymm9
        vpaddw (%rsi, %r13),  %ymm1,  %ymm1
        vpsubw (%rsi, %r14),  %ymm2, %ymm10
        vpaddw (%rsi, %r14),  %ymm2,  %ymm2
        vpsubw (%rsi, %r15),  %ymm3, %ymm11
        vpaddw (%rsi, %r15),  %ymm3,  %ymm3

        mov       12*8( %r9), %r12
        mov        2*8( %r9), %r13
        mov        6*8( %r9), %r14
        mov        1*8( %r9), %r15

        vpsubw (%rsi, %r12),  %ymm4, %ymm12
        vpaddw (%rsi, %r12),  %ymm4,  %ymm4
        vpsubw (%rsi, %r13),  %ymm5, %ymm13
        vpaddw (%rsi, %r13),  %ymm5,  %ymm5
        vpsubw (%rsi, %r14),  %ymm6, %ymm14
        vpaddw (%rsi, %r14),  %ymm6,  %ymm6
        vpsubw (%rsi, %r15),  %ymm7, %ymm15
        vpaddw (%rsi, %r15),  %ymm7,  %ymm7

        vmovdqu       %ymm8, 0x100(%rsp)
        vmovdqu       %ymm9, 0x120(%rsp)
        vmovdqu      %ymm10, 0x140(%rsp)
        vmovdqu      %ymm11, 0x160(%rsp)
        vmovdqu      %ymm12, 0x180(%rsp)
        vmovdqu      %ymm13, 0x1a0(%rsp)
        vmovdqu      %ymm14, 0x1c0(%rsp)
        vmovdqu      %ymm15, 0x1e0(%rsp)

        vpaddw        %ymm4,  %ymm0,  %ymm8
        vpsubw        %ymm4,  %ymm0, %ymm12
        vpaddw        %ymm5,  %ymm1,  %ymm9
        vpsubw        %ymm5,  %ymm1, %ymm13
        vpaddw        %ymm6,  %ymm2, %ymm10
        vpsubw        %ymm6,  %ymm2, %ymm14
        vpaddw        %ymm7,  %ymm3, %ymm11
        vpsubw        %ymm7,  %ymm3, %ymm15

        vmovdqu      %ymm12, 0x080(%rsp)
        vmovdqu      %ymm13, 0x0a0(%rsp)
        vmovdqu      %ymm14, 0x0c0(%rsp)
        vmovdqu      %ymm15, 0x0e0(%rsp)

        vmovdqu 0x000( %r8), %ymm12
        vmovdqu 0x020( %r8), %ymm13

// ================================
// x^4 - 1

        vmovdqu  0x040( %r8),  %ymm13
        barrett_reducex2 0, 1,  8,  9, 12, 13, 14, 15
        barrett_reducex2 2, 3, 10, 11, 12, 13, 14, 15

        vpsubw        %ymm2,  %ymm0, %ymm10
        vpaddw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm1, %ymm11
        vpaddw        %ymm3,  %ymm1,  %ymm1
        vpaddw        %ymm1,  %ymm0,  %ymm8
        vpsubw        %ymm1,  %ymm0,  %ymm9

        vmovdqu       %ymm8, 0x620(%rsp)

        vmovdqu 0x000(%rdx),  %ymm4
        vmovdqu 0x020(%rdx),  %ymm5
        vmovdqu 0x040(%rdx),  %ymm6
        vmovdqu 0x060(%rdx),  %ymm7

        vpmullw 0x000(%rcx),  %ymm8, %ymm14
        vpmullw 0x020(%rcx),  %ymm9, %ymm15
        vpmulhw       %ymm4,  %ymm8,  %ymm4
        vpmulhw       %ymm5,  %ymm9,  %ymm5
        vpmulhw      %ymm14, %ymm12, %ymm14
        vpmulhw      %ymm15, %ymm12, %ymm15
        vpsubw       %ymm14,  %ymm4,  %ymm4
        vpsubw       %ymm15,  %ymm5,  %ymm5

        vpmullw 0x040(%rcx), %ymm10, %ymm14
        vpmullw 0x060(%rcx), %ymm10, %ymm15
        vpmullw 0x040(%rcx), %ymm11,  %ymm8
        vpmullw 0x060(%rcx), %ymm11,  %ymm9
        vpmulhw       %ymm6, %ymm10,  %ymm0
        vpmulhw       %ymm7, %ymm10,  %ymm1
        vpmulhw       %ymm6, %ymm11,  %ymm2
        vpmulhw       %ymm7, %ymm11,  %ymm3
        vpmulhw      %ymm14, %ymm12, %ymm14
        vpmulhw      %ymm15, %ymm12, %ymm15
        vpmulhw       %ymm8, %ymm12,  %ymm8
        vpmulhw       %ymm9, %ymm12,  %ymm9
        vpsubw       %ymm14,  %ymm0,  %ymm0
        vpsubw       %ymm15,  %ymm1,  %ymm1
        vpsubw        %ymm8,  %ymm2,  %ymm2
        vpsubw        %ymm9,  %ymm3,  %ymm3

        vpsubw        %ymm3,  %ymm0, %ymm6
        vpaddw        %ymm2,  %ymm1, %ymm7

        vpaddw        %ymm5,  %ymm4,  %ymm0
        vpsubw        %ymm5,  %ymm4,  %ymm1
        vpsubw        %ymm6,  %ymm0,  %ymm2
        vpaddw        %ymm6,  %ymm0,  %ymm0
        vpsubw        %ymm7,  %ymm1,  %ymm3
        vpaddw        %ymm7,  %ymm1,  %ymm1

        vmovdqu  0x040( %r8),  %ymm13

        barrett_reducex2 0, 1, 0, 1, 12, 13, 14, 15
        barrett_reducex2 2, 3, 2, 3, 12, 13, 14, 15

        vmovdqu       %ymm0, 0x000(%rsp)
        vmovdqu       %ymm1, 0x020(%rsp)
        vmovdqu       %ymm2, 0x040(%rsp)
        vmovdqu       %ymm3, 0x060(%rsp)

        vmovdqu  0x020( %r8),  %ymm13

// ================================
// x^4 + 1

        vmovdqu 0x080(%rsp),  %ymm8
        vmovdqu 0x0a0(%rsp),  %ymm9
        vmovdqu 0x0c0(%rsp), %ymm10
        vmovdqu 0x0e0(%rsp), %ymm11

        vpmullw 0x080(%rcx),  %ymm8,  %ymm0
        vpmullw 0x0a0(%rcx), %ymm11,  %ymm1
        vpmullw 0x0c0(%rcx), %ymm10,  %ymm2
        vpmullw 0x0e0(%rcx),  %ymm9,  %ymm3
        vpmulhw       %ymm0, %ymm12,  %ymm0
        vpmulhw       %ymm1, %ymm12,  %ymm1
        vpmulhw       %ymm2, %ymm12,  %ymm2
        vpmulhw       %ymm3, %ymm12,  %ymm3
        vpmulhw 0x080(%rdx),  %ymm8, %ymm12
        vpmulhw 0x0a0(%rdx), %ymm11, %ymm13
        vpmulhw 0x0c0(%rdx), %ymm10, %ymm14
        vpmulhw 0x0e0(%rdx),  %ymm9, %ymm15
        vpsubw        %ymm0, %ymm12,  %ymm0
        vpsubw        %ymm1, %ymm13,  %ymm1
        vpsubw        %ymm2, %ymm14,  %ymm2
        vpsubw        %ymm3, %ymm15,  %ymm3

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpsubw        %ymm1,  %ymm0,  %ymm0
        vpsubw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm0,  %ymm0

        barrett_reduce 0, 0, 12, 13, 14
        vmovdqu       %ymm0, 0x080(%rsp)

        vpmullw 0x080(%rcx),  %ymm9,  %ymm4
        vpmullw 0x0a0(%rcx),  %ymm8,  %ymm5
        vpmullw 0x0c0(%rcx), %ymm11,  %ymm6
        vpmullw 0x0e0(%rcx), %ymm10,  %ymm7
        vpmulhw       %ymm4, %ymm12,  %ymm4
        vpmulhw       %ymm5, %ymm12,  %ymm5
        vpmulhw       %ymm6, %ymm12,  %ymm6
        vpmulhw       %ymm7, %ymm12,  %ymm7
        vpmulhw 0x080(%rdx),  %ymm9, %ymm12
        vpmulhw 0x0a0(%rdx),  %ymm8, %ymm13
        vpmulhw 0x0c0(%rdx), %ymm11, %ymm14
        vpmulhw 0x0e0(%rdx), %ymm10, %ymm15
        vpsubw        %ymm4, %ymm12,  %ymm4
        vpsubw        %ymm5, %ymm13,  %ymm5
        vpsubw        %ymm6, %ymm14,  %ymm6
        vpsubw        %ymm7, %ymm15,  %ymm7

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm5,  %ymm4,  %ymm4
        vpsubw        %ymm6,  %ymm4,  %ymm4
        vpsubw        %ymm7,  %ymm4,  %ymm4

        barrett_reduce 4, 4, 12, 13, 14
        vmovdqu       %ymm4, 0x0a0(%rsp)

        vpmullw 0x080(%rcx), %ymm10,  %ymm0
        vpmullw 0x0a0(%rcx),  %ymm9,  %ymm1
        vpmullw 0x0c0(%rcx),  %ymm8,  %ymm2
        vpmullw 0x0e0(%rcx), %ymm11,  %ymm3
        vpmulhw       %ymm0, %ymm12,  %ymm0
        vpmulhw       %ymm1, %ymm12,  %ymm1
        vpmulhw       %ymm2, %ymm12,  %ymm2
        vpmulhw       %ymm3, %ymm12,  %ymm3
        vpmulhw 0x080(%rdx), %ymm10, %ymm12
        vpmulhw 0x0a0(%rdx),  %ymm9, %ymm13
        vpmulhw 0x0c0(%rdx),  %ymm8, %ymm14
        vpmulhw 0x0e0(%rdx), %ymm11, %ymm15
        vpsubw        %ymm0, %ymm12,  %ymm0
        vpsubw        %ymm1, %ymm13,  %ymm1
        vpsubw        %ymm2, %ymm14,  %ymm2
        vpsubw        %ymm3, %ymm15,  %ymm3

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm1,  %ymm0,  %ymm0
        vpaddw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm0,  %ymm0

        barrett_reduce 0, 0, 12, 13, 14
        vmovdqu       %ymm0, 0x0c0(%rsp)

        vpmullw 0x080(%rcx), %ymm11,  %ymm4
        vpmullw 0x0a0(%rcx), %ymm10,  %ymm5
        vpmullw 0x0c0(%rcx),  %ymm9,  %ymm6
        vpmullw 0x0e0(%rcx),  %ymm8,  %ymm7
        vpmulhw       %ymm4, %ymm12,  %ymm4
        vpmulhw       %ymm5, %ymm12,  %ymm5
        vpmulhw       %ymm6, %ymm12,  %ymm6
        vpmulhw       %ymm7, %ymm12,  %ymm7
        vpmulhw 0x080(%rdx), %ymm11, %ymm12
        vpmulhw 0x0a0(%rdx), %ymm10, %ymm13
        vpmulhw 0x0c0(%rdx),  %ymm9, %ymm14
        vpmulhw 0x0e0(%rdx),  %ymm8, %ymm15
        vpsubw        %ymm4, %ymm12,  %ymm4
        vpsubw        %ymm5, %ymm13,  %ymm5
        vpsubw        %ymm6, %ymm14,  %ymm6
        vpsubw        %ymm7, %ymm15,  %ymm7

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm5,  %ymm4,  %ymm4
        vpaddw        %ymm6,  %ymm4,  %ymm4
        vpaddw        %ymm7,  %ymm4,  %ymm4

        barrett_reduce 4, 4, 12, 13, 14
        vmovdqu       %ymm4, 0x0e0(%rsp)

// ================================
// x^8 + 1

        __asm_negacyclic_karatsuba8_precompute_core \
                %rsp, %rsp, %rdx, %rcx, %rsp, 0x200, 0x100, 0x100, 0x100, 0x300, \
                0x000( %r8), 0x020( %r8), 0x040( %r8)

        vmovdqu 0x000(%rsp),  %ymm0
        vmovdqu 0x020(%rsp),  %ymm1
        vmovdqu 0x040(%rsp),  %ymm2
        vmovdqu 0x060(%rsp),  %ymm3

        vpsubw  0x080(%rsp),  %ymm0,  %ymm4
        vpaddw  0x080(%rsp),  %ymm0,  %ymm0
        vpsubw  0x0a0(%rsp),  %ymm1,  %ymm5
        vpaddw  0x0a0(%rsp),  %ymm1,  %ymm1
        vpsubw  0x0c0(%rsp),  %ymm2,  %ymm6
        vpaddw  0x0c0(%rsp),  %ymm2,  %ymm2
        vpsubw  0x0e0(%rsp),  %ymm3,  %ymm7
        vpaddw  0x0e0(%rsp),  %ymm3,  %ymm3

        vpsubw  0x200(%rsp),  %ymm0,  %ymm8
        vpaddw  0x200(%rsp),  %ymm0,  %ymm0
        vpsubw  0x220(%rsp),  %ymm1,  %ymm9
        vpaddw  0x220(%rsp),  %ymm1,  %ymm1
        vpsubw  0x240(%rsp),  %ymm2, %ymm10
        vpaddw  0x240(%rsp),  %ymm2,  %ymm2
        vpsubw  0x260(%rsp),  %ymm3, %ymm11
        vpaddw  0x260(%rsp),  %ymm3,  %ymm3
        vpsubw  0x280(%rsp),  %ymm4, %ymm12
        vpaddw  0x280(%rsp),  %ymm4,  %ymm4
        vpsubw  0x2a0(%rsp),  %ymm5, %ymm13
        vpaddw  0x2a0(%rsp),  %ymm5,  %ymm5
        vpsubw  0x2c0(%rsp),  %ymm6, %ymm14
        vpaddw  0x2c0(%rsp),  %ymm6,  %ymm6
        vpsubw  0x2e0(%rsp),  %ymm7, %ymm15
        vpaddw  0x2e0(%rsp),  %ymm7,  %ymm7

        vpaddw  0x600(%rsp),  %ymm0,  %ymm0
        vpaddw  0x600(%rsp),  %ymm1,  %ymm1
        vpaddw  0x600(%rsp),  %ymm2,  %ymm2
        vpaddw  0x600(%rsp),  %ymm3,  %ymm3
        vpaddw  0x600(%rsp),  %ymm4,  %ymm4
        vpaddw  0x600(%rsp),  %ymm5,  %ymm5
        vpaddw  0x600(%rsp),  %ymm6,  %ymm6
        vpaddw  0x600(%rsp),  %ymm7,  %ymm7
        vpaddw  0x600(%rsp),  %ymm8,  %ymm8
        vpaddw  0x600(%rsp),  %ymm9,  %ymm9
        vpaddw  0x600(%rsp), %ymm10, %ymm10
        vpaddw  0x600(%rsp), %ymm11, %ymm11
        vpaddw  0x600(%rsp), %ymm12, %ymm12
        vpaddw  0x600(%rsp), %ymm13, %ymm13
        vpaddw  0x600(%rsp), %ymm14, %ymm14
        vpaddw  0x600(%rsp), %ymm15, %ymm15


.if 0

// g = 3
// map x to g^{16 - x}
size_t rader_out_permute[17] = {
0,
6, 2, 12, 4,
7, 8, 14, 16,
11, 15, 5, 13,
10, 9, 3, 1
};

.endif

// 0, 576, 192, 1152, 384, 672, 768, 1344, 1536, 1056, 1440, 480, 1248, 960, 864, 288, 96

        vmovdqu       %ymm0,   576(%rdi)
        vmovdqu       %ymm1,   192(%rdi)
        vmovdqu       %ymm2,  1152(%rdi)
        vmovdqu       %ymm3,   384(%rdi)
        vmovdqu       %ymm4,   672(%rdi)
        vmovdqu       %ymm5,   768(%rdi)
        vmovdqu       %ymm6,  1344(%rdi)
        vmovdqu       %ymm7,  1536(%rdi)
        vmovdqu       %ymm8,  1056(%rdi)
        vmovdqu       %ymm9,  1440(%rdi)
        vmovdqu      %ymm10,   480(%rdi)
        vmovdqu      %ymm11,  1248(%rdi)
        vmovdqu      %ymm12,   960(%rdi)
        vmovdqu      %ymm13,   864(%rdi)
        vmovdqu      %ymm14,   288(%rdi)
        vmovdqu      %ymm15,    96(%rdi)

        vmovdqu 0x620(%rsp),  %ymm0
        vpaddw  0x600(%rsp),  %ymm0,  %ymm0
        vmovdqu       %ymm0,      0(%rdi)


        vzeroupper
        mov    %rbp, %rsp
        pop    %r15
        pop    %r14
        pop    %r13
        pop    %r12
        pop    %rbp
        ret

.global __asm_rader17_post_tbl
.global ___asm_rader17_post_tbl
__asm_rader17_post_tbl:
___asm_rader17_post_tbl:

        push   %rbp
        push   %r12
        push   %r13
        push   %r14
        push   %r15
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x700, %rsp

        vmovdqu     0(%rsi),  %ymm0
        vmovdqu       %ymm0, 0x600(%rsp)
        add     $0x020, %rdx
        add     $0x020, %rcx

.if 0

// g = 3
// map x to g^x
size_t rader_in_permute[17] = {
0,
3, 9, 10, 13,
5, 15, 11, 16,
14, 8, 7, 4,
12, 2, 6, 1
};

.endif

// 0, 288, 864, 960, 1248, 480, 1440, 1056, 1536, 1344, 768, 672, 384, 1152, 192, 576, 96

        vmovdqu   288(%rsi),  %ymm0
        vmovdqu   864(%rsi),  %ymm1
        vmovdqu   960(%rsi),  %ymm2
        vmovdqu  1248(%rsi),  %ymm3
        vmovdqu   480(%rsi),  %ymm4
        vmovdqu  1440(%rsi),  %ymm5
        vmovdqu  1056(%rsi),  %ymm6
        vmovdqu  1536(%rsi),  %ymm7

        vpsubw   1344(%rsi),  %ymm0,  %ymm8
        vpaddw   1344(%rsi),  %ymm0,  %ymm0
        vpsubw    768(%rsi),  %ymm1,  %ymm9
        vpaddw    768(%rsi),  %ymm1,  %ymm1
        vpsubw    672(%rsi),  %ymm2, %ymm10
        vpaddw    672(%rsi),  %ymm2,  %ymm2
        vpsubw    384(%rsi),  %ymm3, %ymm11
        vpaddw    384(%rsi),  %ymm3,  %ymm3
        vpsubw   1152(%rsi),  %ymm4, %ymm12
        vpaddw   1152(%rsi),  %ymm4,  %ymm4
        vpsubw    192(%rsi),  %ymm5, %ymm13
        vpaddw    192(%rsi),  %ymm5,  %ymm5
        vpsubw    576(%rsi),  %ymm6, %ymm14
        vpaddw    576(%rsi),  %ymm6,  %ymm6
        vpsubw     96(%rsi),  %ymm7, %ymm15
        vpaddw     96(%rsi),  %ymm7,  %ymm7

        vmovdqu       %ymm8, 0x100(%rsp)
        vmovdqu       %ymm9, 0x120(%rsp)
        vmovdqu      %ymm10, 0x140(%rsp)
        vmovdqu      %ymm11, 0x160(%rsp)
        vmovdqu      %ymm12, 0x180(%rsp)
        vmovdqu      %ymm13, 0x1a0(%rsp)
        vmovdqu      %ymm14, 0x1c0(%rsp)
        vmovdqu      %ymm15, 0x1e0(%rsp)

        vpaddw        %ymm4,  %ymm0,  %ymm8
        vpsubw        %ymm4,  %ymm0, %ymm12
        vpaddw        %ymm5,  %ymm1,  %ymm9
        vpsubw        %ymm5,  %ymm1, %ymm13
        vpaddw        %ymm6,  %ymm2, %ymm10
        vpsubw        %ymm6,  %ymm2, %ymm14
        vpaddw        %ymm7,  %ymm3, %ymm11
        vpsubw        %ymm7,  %ymm3, %ymm15

        vmovdqu      %ymm12, 0x080(%rsp)
        vmovdqu      %ymm13, 0x0a0(%rsp)
        vmovdqu      %ymm14, 0x0c0(%rsp)
        vmovdqu      %ymm15, 0x0e0(%rsp)

        vmovdqu 0x000( %r8), %ymm12
        vmovdqu 0x020( %r8), %ymm13

// ================================
// x^4 - 1

        vmovdqu  0x040( %r8),  %ymm13
        barrett_reducex2 0, 1,  8,  9, 12, 13, 14, 15
        barrett_reducex2 2, 3, 10, 11, 12, 13, 14, 15

        vpsubw        %ymm2,  %ymm0, %ymm10
        vpaddw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm1, %ymm11
        vpaddw        %ymm3,  %ymm1,  %ymm1
        vpaddw        %ymm1,  %ymm0,  %ymm8
        vpsubw        %ymm1,  %ymm0,  %ymm9

        vmovdqu       %ymm8, 0x620(%rsp)

        vmovdqu 0x000(%rdx),  %ymm4
        vmovdqu 0x020(%rdx),  %ymm5
        vmovdqu 0x040(%rdx),  %ymm6
        vmovdqu 0x060(%rdx),  %ymm7

        vpmullw 0x000(%rcx),  %ymm8, %ymm14
        vpmullw 0x020(%rcx),  %ymm9, %ymm15
        vpmulhw       %ymm4,  %ymm8,  %ymm4
        vpmulhw       %ymm5,  %ymm9,  %ymm5
        vpmulhw      %ymm14, %ymm12, %ymm14
        vpmulhw      %ymm15, %ymm12, %ymm15
        vpsubw       %ymm14,  %ymm4,  %ymm4
        vpsubw       %ymm15,  %ymm5,  %ymm5

        vpmullw 0x040(%rcx), %ymm10, %ymm14
        vpmullw 0x060(%rcx), %ymm10, %ymm15
        vpmullw 0x040(%rcx), %ymm11,  %ymm8
        vpmullw 0x060(%rcx), %ymm11,  %ymm9
        vpmulhw       %ymm6, %ymm10,  %ymm0
        vpmulhw       %ymm7, %ymm10,  %ymm1
        vpmulhw       %ymm6, %ymm11,  %ymm2
        vpmulhw       %ymm7, %ymm11,  %ymm3
        vpmulhw      %ymm14, %ymm12, %ymm14
        vpmulhw      %ymm15, %ymm12, %ymm15
        vpmulhw       %ymm8, %ymm12,  %ymm8
        vpmulhw       %ymm9, %ymm12,  %ymm9
        vpsubw       %ymm14,  %ymm0,  %ymm0
        vpsubw       %ymm15,  %ymm1,  %ymm1
        vpsubw        %ymm8,  %ymm2,  %ymm2
        vpsubw        %ymm9,  %ymm3,  %ymm3

        vpsubw        %ymm3,  %ymm0, %ymm6
        vpaddw        %ymm2,  %ymm1, %ymm7

        vpaddw        %ymm5,  %ymm4,  %ymm0
        vpsubw        %ymm5,  %ymm4,  %ymm1
        vpsubw        %ymm6,  %ymm0,  %ymm2
        vpaddw        %ymm6,  %ymm0,  %ymm0
        vpsubw        %ymm7,  %ymm1,  %ymm3
        vpaddw        %ymm7,  %ymm1,  %ymm1

        vmovdqu  0x040( %r8),  %ymm13

        barrett_reducex2 0, 1, 0, 1, 12, 13, 14, 15
        barrett_reducex2 2, 3, 2, 3, 12, 13, 14, 15

        vmovdqu       %ymm0, 0x000(%rsp)
        vmovdqu       %ymm1, 0x020(%rsp)
        vmovdqu       %ymm2, 0x040(%rsp)
        vmovdqu       %ymm3, 0x060(%rsp)

        vmovdqu  0x020( %r8),  %ymm13

// ================================
// x^4 + 1

        vmovdqu 0x080(%rsp),  %ymm8
        vmovdqu 0x0a0(%rsp),  %ymm9
        vmovdqu 0x0c0(%rsp), %ymm10
        vmovdqu 0x0e0(%rsp), %ymm11

        vpmullw 0x080(%rcx),  %ymm8,  %ymm0
        vpmullw 0x0a0(%rcx), %ymm11,  %ymm1
        vpmullw 0x0c0(%rcx), %ymm10,  %ymm2
        vpmullw 0x0e0(%rcx),  %ymm9,  %ymm3
        vpmulhw       %ymm0, %ymm12,  %ymm0
        vpmulhw       %ymm1, %ymm12,  %ymm1
        vpmulhw       %ymm2, %ymm12,  %ymm2
        vpmulhw       %ymm3, %ymm12,  %ymm3
        vpmulhw 0x080(%rdx),  %ymm8, %ymm12
        vpmulhw 0x0a0(%rdx), %ymm11, %ymm13
        vpmulhw 0x0c0(%rdx), %ymm10, %ymm14
        vpmulhw 0x0e0(%rdx),  %ymm9, %ymm15
        vpsubw        %ymm0, %ymm12,  %ymm0
        vpsubw        %ymm1, %ymm13,  %ymm1
        vpsubw        %ymm2, %ymm14,  %ymm2
        vpsubw        %ymm3, %ymm15,  %ymm3

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpsubw        %ymm1,  %ymm0,  %ymm0
        vpsubw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm0,  %ymm0

        barrett_reduce 0, 0, 12, 13, 14
        vmovdqu       %ymm0, 0x080(%rsp)

        vpmullw 0x080(%rcx),  %ymm9,  %ymm4
        vpmullw 0x0a0(%rcx),  %ymm8,  %ymm5
        vpmullw 0x0c0(%rcx), %ymm11,  %ymm6
        vpmullw 0x0e0(%rcx), %ymm10,  %ymm7
        vpmulhw       %ymm4, %ymm12,  %ymm4
        vpmulhw       %ymm5, %ymm12,  %ymm5
        vpmulhw       %ymm6, %ymm12,  %ymm6
        vpmulhw       %ymm7, %ymm12,  %ymm7
        vpmulhw 0x080(%rdx),  %ymm9, %ymm12
        vpmulhw 0x0a0(%rdx),  %ymm8, %ymm13
        vpmulhw 0x0c0(%rdx), %ymm11, %ymm14
        vpmulhw 0x0e0(%rdx), %ymm10, %ymm15
        vpsubw        %ymm4, %ymm12,  %ymm4
        vpsubw        %ymm5, %ymm13,  %ymm5
        vpsubw        %ymm6, %ymm14,  %ymm6
        vpsubw        %ymm7, %ymm15,  %ymm7

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm5,  %ymm4,  %ymm4
        vpsubw        %ymm6,  %ymm4,  %ymm4
        vpsubw        %ymm7,  %ymm4,  %ymm4

        barrett_reduce 4, 4, 12, 13, 14
        vmovdqu       %ymm4, 0x0a0(%rsp)

        vpmullw 0x080(%rcx), %ymm10,  %ymm0
        vpmullw 0x0a0(%rcx),  %ymm9,  %ymm1
        vpmullw 0x0c0(%rcx),  %ymm8,  %ymm2
        vpmullw 0x0e0(%rcx), %ymm11,  %ymm3
        vpmulhw       %ymm0, %ymm12,  %ymm0
        vpmulhw       %ymm1, %ymm12,  %ymm1
        vpmulhw       %ymm2, %ymm12,  %ymm2
        vpmulhw       %ymm3, %ymm12,  %ymm3
        vpmulhw 0x080(%rdx), %ymm10, %ymm12
        vpmulhw 0x0a0(%rdx),  %ymm9, %ymm13
        vpmulhw 0x0c0(%rdx),  %ymm8, %ymm14
        vpmulhw 0x0e0(%rdx), %ymm11, %ymm15
        vpsubw        %ymm0, %ymm12,  %ymm0
        vpsubw        %ymm1, %ymm13,  %ymm1
        vpsubw        %ymm2, %ymm14,  %ymm2
        vpsubw        %ymm3, %ymm15,  %ymm3

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm1,  %ymm0,  %ymm0
        vpaddw        %ymm2,  %ymm0,  %ymm0
        vpsubw        %ymm3,  %ymm0,  %ymm0

        barrett_reduce 0, 0, 12, 13, 14
        vmovdqu       %ymm0, 0x0c0(%rsp)

        vpmullw 0x080(%rcx), %ymm11,  %ymm4
        vpmullw 0x0a0(%rcx), %ymm10,  %ymm5
        vpmullw 0x0c0(%rcx),  %ymm9,  %ymm6
        vpmullw 0x0e0(%rcx),  %ymm8,  %ymm7
        vpmulhw       %ymm4, %ymm12,  %ymm4
        vpmulhw       %ymm5, %ymm12,  %ymm5
        vpmulhw       %ymm6, %ymm12,  %ymm6
        vpmulhw       %ymm7, %ymm12,  %ymm7
        vpmulhw 0x080(%rdx), %ymm11, %ymm12
        vpmulhw 0x0a0(%rdx), %ymm10, %ymm13
        vpmulhw 0x0c0(%rdx),  %ymm9, %ymm14
        vpmulhw 0x0e0(%rdx),  %ymm8, %ymm15
        vpsubw        %ymm4, %ymm12,  %ymm4
        vpsubw        %ymm5, %ymm13,  %ymm5
        vpsubw        %ymm6, %ymm14,  %ymm6
        vpsubw        %ymm7, %ymm15,  %ymm7

        vmovdqu  0x000( %r8),  %ymm12
        vmovdqu  0x040( %r8),  %ymm13

        vpaddw        %ymm5,  %ymm4,  %ymm4
        vpaddw        %ymm6,  %ymm4,  %ymm4
        vpaddw        %ymm7,  %ymm4,  %ymm4

        barrett_reduce 4, 4, 12, 13, 14
        vmovdqu       %ymm4, 0x0e0(%rsp)

// ================================
// x^8 + 1

        __asm_negacyclic_karatsuba8_precompute_core \
                %rsp, %rsp, %rdx, %rcx, %rsp, 0x200, 0x100, 0x100, 0x100, 0x300, \
                0x000( %r8), 0x020( %r8), 0x040( %r8)

        vmovdqu 0x000(%rsp),  %ymm0
        vmovdqu 0x020(%rsp),  %ymm1
        vmovdqu 0x040(%rsp),  %ymm2
        vmovdqu 0x060(%rsp),  %ymm3

        vpsubw  0x080(%rsp),  %ymm0,  %ymm4
        vpaddw  0x080(%rsp),  %ymm0,  %ymm0
        vpsubw  0x0a0(%rsp),  %ymm1,  %ymm5
        vpaddw  0x0a0(%rsp),  %ymm1,  %ymm1
        vpsubw  0x0c0(%rsp),  %ymm2,  %ymm6
        vpaddw  0x0c0(%rsp),  %ymm2,  %ymm2
        vpsubw  0x0e0(%rsp),  %ymm3,  %ymm7
        vpaddw  0x0e0(%rsp),  %ymm3,  %ymm3

        vpsubw  0x200(%rsp),  %ymm0,  %ymm8
        vpaddw  0x200(%rsp),  %ymm0,  %ymm0
        vpsubw  0x220(%rsp),  %ymm1,  %ymm9
        vpaddw  0x220(%rsp),  %ymm1,  %ymm1
        vpsubw  0x240(%rsp),  %ymm2, %ymm10
        vpaddw  0x240(%rsp),  %ymm2,  %ymm2
        vpsubw  0x260(%rsp),  %ymm3, %ymm11
        vpaddw  0x260(%rsp),  %ymm3,  %ymm3
        vpsubw  0x280(%rsp),  %ymm4, %ymm12
        vpaddw  0x280(%rsp),  %ymm4,  %ymm4
        vpsubw  0x2a0(%rsp),  %ymm5, %ymm13
        vpaddw  0x2a0(%rsp),  %ymm5,  %ymm5
        vpsubw  0x2c0(%rsp),  %ymm6, %ymm14
        vpaddw  0x2c0(%rsp),  %ymm6,  %ymm6
        vpsubw  0x2e0(%rsp),  %ymm7, %ymm15
        vpaddw  0x2e0(%rsp),  %ymm7,  %ymm7

        vpaddw  0x600(%rsp),  %ymm0,  %ymm0
        vpaddw  0x600(%rsp),  %ymm1,  %ymm1
        vpaddw  0x600(%rsp),  %ymm2,  %ymm2
        vpaddw  0x600(%rsp),  %ymm3,  %ymm3
        vpaddw  0x600(%rsp),  %ymm4,  %ymm4
        vpaddw  0x600(%rsp),  %ymm5,  %ymm5
        vpaddw  0x600(%rsp),  %ymm6,  %ymm6
        vpaddw  0x600(%rsp),  %ymm7,  %ymm7
        vpaddw  0x600(%rsp),  %ymm8,  %ymm8
        vpaddw  0x600(%rsp),  %ymm9,  %ymm9
        vpaddw  0x600(%rsp), %ymm10, %ymm10
        vpaddw  0x600(%rsp), %ymm11, %ymm11
        vpaddw  0x600(%rsp), %ymm12, %ymm12
        vpaddw  0x600(%rsp), %ymm13, %ymm13
        vpaddw  0x600(%rsp), %ymm14, %ymm14
        vpaddw  0x600(%rsp), %ymm15, %ymm15


.if 0

// g = 3
// map x to g^{16 - x}
size_t rader_out_permute[17] = {
0,
6, 2, 12, 4,
7, 8, 14, 16,
11, 15, 5, 13,
10, 9, 3, 1
};

.endif

// 0, 576, 192, 1152, 384, 672, 768, 1344, 1536, 1056, 1440, 480, 1248, 960, 864, 288, 96

        mov       6*8( %r9), %r12
        mov       2*8( %r9), %r13
        mov      12*8( %r9), %r14
        mov       4*8( %r9), %r15

        vmovdqu       %ymm0, (%rdi, %r12)
        vmovdqu       %ymm1, (%rdi, %r13)
        vmovdqu       %ymm2, (%rdi, %r14)
        vmovdqu       %ymm3, (%rdi, %r15)

        mov       7*8( %r9), %r12
        mov       8*8( %r9), %r13
        mov      14*8( %r9), %r14
        mov      16*8( %r9), %r15

        vmovdqu       %ymm4, (%rdi, %r12)
        vmovdqu       %ymm5, (%rdi, %r13)
        vmovdqu       %ymm6, (%rdi, %r14)
        vmovdqu       %ymm7, (%rdi, %r15)

        mov      11*8( %r9), %r12
        mov      15*8( %r9), %r13
        mov       5*8( %r9), %r14
        mov      13*8( %r9), %r15

        vmovdqu       %ymm8, (%rdi, %r12)
        vmovdqu       %ymm9, (%rdi, %r13)
        vmovdqu      %ymm10, (%rdi, %r14)
        vmovdqu      %ymm11, (%rdi, %r15)

        mov      10*8( %r9), %r12
        mov       9*8( %r9), %r13
        mov       3*8( %r9), %r14
        mov       1*8( %r9), %r15

        vmovdqu      %ymm12, (%rdi, %r12)
        vmovdqu      %ymm13, (%rdi, %r13)
        vmovdqu      %ymm14, (%rdi, %r14)
        vmovdqu      %ymm15, (%rdi, %r15)

        vmovdqu 0x620(%rsp),  %ymm0
        vpaddw  0x600(%rsp),  %ymm0,  %ymm0

        mov       0*8( %r9), %r12

        vmovdqu       %ymm0, (%rdi, %r12)


        vzeroupper
        mov    %rbp, %rsp
        pop    %r15
        pop    %r14
        pop    %r13
        pop    %r12
        pop    %rbp
        ret





