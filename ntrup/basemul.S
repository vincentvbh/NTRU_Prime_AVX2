

#include "basemul_core.inc"

#include "permute.inc"

.text

.global __asm_permute_test
.global ___asm_permute_test
__asm_permute_test:
___asm_permute_test:

        vmovdqu 0x00(%rsi), %ymm0
        vmovdqu 0x20(%rsi), %ymm1

        __asm_interleave 2, 3, 0, 1

        vpxor %ymm8, %ymm8, %ymm8

        __asm_deinterleave 12, 13, 2, 3, 4, 5, 8

        vmovdqu      %ymm12, 0x00(%rdi)
        vmovdqu      %ymm13, 0x20(%rdi)

        ret

.global __asm_schoolbook2
.global ___asm_schoolbook2
__asm_schoolbook2:
___asm_schoolbook2:

        vmovdqu     (%rsi),  %ymm0
        vmovdqu 0x20(%rsi),  %ymm1

        vmovdqu     (%rdx),  %ymm2
        vmovdqu 0x20(%rdx),  %ymm3

        vmovdqu     (%rcx), %ymm12
        vmovdqu 0x20(%rcx), %ymm13

        __asm_schoolbook2_core %rdi, 0x00, \
                0, 1, 2, 3, \
                12, 13, \
                8, 9, 10, 11, 14, 15

        ret

.global __asm_cyclic_schoolbook2
.global ___asm_cyclic_schoolbook2
__asm_cyclic_schoolbook2:
___asm_cyclic_schoolbook2:

        vmovdqu     (%rsi),  %ymm0
        vmovdqu 0x20(%rsi),  %ymm1

        vmovdqu     (%rdx),  %ymm2
        vmovdqu 0x20(%rdx),  %ymm3

        vmovdqu     (%rcx), %ymm12
        vmovdqu 0x20(%rcx), %ymm13

        __asm_cyclic_schoolbook2_core %rdi, 0x00, \
                0, 1, 2, 3, \
                12, 13, \
                8, 9, 10, 11, 14, 15

        ret

.global __asm_negacyclic_schoolbook2
.global ___asm_negacyclic_schoolbook2
__asm_negacyclic_schoolbook2:
___asm_negacyclic_schoolbook2:

        vmovdqu     (%rsi),  %ymm0
        vmovdqu 0x20(%rsi),  %ymm1

        vmovdqu     (%rdx),  %ymm2
        vmovdqu 0x20(%rdx),  %ymm3

        vmovdqu     (%rcx), %ymm12
        vmovdqu 0x20(%rcx), %ymm13

        __asm_negacyclic_schoolbook2_core %rdi, 0x00, \
                0, 1, 2, 3, \
                12, 13, \
                8, 9, 10, 11, 14, 15

        ret


.global __asm_schoolbook4
.global ___asm_schoolbook4
__asm_schoolbook4:
___asm_schoolbook4:

        vmovdqu     (%rsi),  %ymm0
        vmovdqu 0x20(%rsi),  %ymm1
        vmovdqu 0x40(%rsi),  %ymm2
        vmovdqu 0x60(%rsi),  %ymm3
        vmovdqu     (%rdx),  %ymm4
        vmovdqu 0x20(%rdx),  %ymm5
        vmovdqu 0x40(%rdx),  %ymm6
        vmovdqu 0x60(%rdx),  %ymm7

        vmovdqu     (%rcx), %ymm12
        vmovdqu 0x20(%rcx), %ymm13

        __asm_schoolbook4_core %rdi, 0x00, \
                0, 1, 2, 3, \
                4, 5, 6, 7, \
                12, 13, 0x20(%rcx), 0x40(%rcx), \
                8, 9, 10, 11, 14, 15

        ret

.global __asm_cyclic_schoolbook4
.global ___asm_cyclic_schoolbook4
__asm_cyclic_schoolbook4:
___asm_cyclic_schoolbook4:

        vmovdqu     (%rsi),  %ymm0
        vmovdqu 0x20(%rsi),  %ymm1
        vmovdqu 0x40(%rsi),  %ymm2
        vmovdqu 0x60(%rsi),  %ymm3

        vmovdqu     (%rdx),  %ymm4
        vmovdqu 0x20(%rdx),  %ymm5
        vmovdqu 0x40(%rdx),  %ymm6
        vmovdqu 0x60(%rdx),  %ymm7

        vmovdqu     (%rcx), %ymm12
        vmovdqu 0x20(%rcx), %ymm13

        __asm_cyclic_schoolbook4_core %rdi, 0x00, \
                0, 1, 2, 3, \
                4, 5, 6, 7, \
                12, 13, 0x20(%rcx), 0x40(%rcx), \
                8, 9, 10, 11, 14, 15

        ret

.global __asm_negacyclic_schoolbook4
.global ___asm_negacyclic_schoolbook4
__asm_negacyclic_schoolbook4:
___asm_negacyclic_schoolbook4:

        vmovdqu     (%rsi),  %ymm0
        vmovdqu 0x20(%rsi),  %ymm1
        vmovdqu 0x40(%rsi),  %ymm2
        vmovdqu 0x60(%rsi),  %ymm3

        vmovdqu     (%rdx),  %ymm4
        vmovdqu 0x20(%rdx),  %ymm5
        vmovdqu 0x40(%rdx),  %ymm6
        vmovdqu 0x60(%rdx),  %ymm7

        vmovdqu     (%rcx), %ymm12
        vmovdqu 0x20(%rcx), %ymm13

        __asm_negacyclic_schoolbook4_core %rdi, 0x00, \
                0, 1, 2, 3, \
                4, 5, 6, 7, \
                12, 13, 0x20(%rcx), 0x40(%rcx), \
                8, 9, 10, 11, 14, 15

        ret




.global __asm_cyclic_FFT2
.global ___asm_cyclic_FFT2
__asm_cyclic_FFT2:
___asm_cyclic_FFT2:

        vmovdqu     (%rsi),  %ymm0
        vmovdqu 0x20(%rsi),  %ymm1
        vmovdqu     (%rdx),  %ymm2
        vmovdqu 0x20(%rdx),  %ymm3

        vmovdqu     (%rcx), %ymm12
        vmovdqu 0x20(%rcx), %ymm13

        vpaddw       %ymm1,  %ymm0,  %ymm4
        vpsubw       %ymm1,  %ymm0,  %ymm5
        vpaddw       %ymm3,  %ymm2,  %ymm6
        vpsubw       %ymm3,  %ymm2,  %ymm7

        vpmullw      %ymm4,  %ymm6,  %ymm0
        vpmulhw      %ymm4,  %ymm6,  %ymm1
        vpmullw      %ymm5,  %ymm7,  %ymm2
        vpmulhw      %ymm5,  %ymm7,  %ymm3
        vpmullw     %ymm13,  %ymm0,  %ymm0
        vpmullw     %ymm13,  %ymm2,  %ymm2
        vpmulhw     %ymm12,  %ymm0,  %ymm0
        vpmulhw     %ymm12,  %ymm2,  %ymm2
        vpsubw       %ymm0,  %ymm1,  %ymm1
        vpsubw       %ymm2,  %ymm3,  %ymm3
        vpaddw       %ymm3,  %ymm1,  %ymm4
        vpsubw       %ymm3,  %ymm1,  %ymm6

        vmovdqu      %ymm4,     (%rdi)
        vmovdqu      %ymm6, 0x20(%rdi)

        ret

.global __asm_cyclic_FFT4
.global ___asm_cyclic_FFT4
__asm_cyclic_FFT4:
___asm_cyclic_FFT4:

        vmovdqu 0x00(%rsi), %ymm0
        vmovdqu 0x20(%rsi), %ymm1
        vmovdqu 0x40(%rsi), %ymm2
        vmovdqu 0x60(%rsi), %ymm3

        vmovdqu 0x00(%rdx), %ymm4
        vmovdqu 0x20(%rdx), %ymm5
        vmovdqu 0x40(%rdx), %ymm6
        vmovdqu 0x60(%rdx), %ymm7


        vpaddw       %ymm2,  %ymm0,  %ymm8
        vpsubw       %ymm2,  %ymm0, %ymm10
        vpaddw       %ymm3,  %ymm1,  %ymm9
        vpsubw       %ymm3,  %ymm1, %ymm11

        vpaddw       %ymm6,  %ymm4, %ymm12
        vpsubw       %ymm6,  %ymm4, %ymm14
        vpaddw       %ymm7,  %ymm5, %ymm13
        vpsubw       %ymm7,  %ymm5, %ymm15

// ================================

        vpmullw     %ymm12,  %ymm8,  %ymm4
        vpmulhw     %ymm12,  %ymm8,  %ymm6
        vpmullw     %ymm13,  %ymm9,  %ymm5
        vpmulhw     %ymm13,  %ymm9,  %ymm7

        vpmullw 0x20(%rcx),  %ymm4,  %ymm4
        vpmullw 0x20(%rcx),  %ymm5,  %ymm5
        vpmulhw 0x00(%rcx),  %ymm4,  %ymm4
        vpmulhw 0x00(%rcx),  %ymm5,  %ymm5

        vpsubw       %ymm4,  %ymm6,  %ymm6
        vpsubw       %ymm5,  %ymm7,  %ymm7

        vpaddw       %ymm7,  %ymm6,  %ymm2

// ================

        vpmullw     %ymm13,  %ymm8,  %ymm4
        vpmulhw     %ymm13,  %ymm8,  %ymm6
        vpmullw     %ymm12,  %ymm9,  %ymm5
        vpmulhw     %ymm12,  %ymm9,  %ymm7

        vpmullw 0x20(%rcx),  %ymm4,  %ymm4
        vpmullw 0x20(%rcx),  %ymm5,  %ymm5
        vpmulhw 0x00(%rcx),  %ymm4,  %ymm4
        vpmulhw 0x00(%rcx),  %ymm5,  %ymm5

        vpsubw       %ymm4,  %ymm6,  %ymm6
        vpsubw       %ymm5,  %ymm7,  %ymm7

        vpaddw       %ymm7,  %ymm6,  %ymm3

// ================================

        vpmullw     %ymm14, %ymm10,  %ymm4
        vpmulhw     %ymm14, %ymm10,  %ymm6
        vpmullw     %ymm15, %ymm11,  %ymm5
        vpmulhw     %ymm15, %ymm11,  %ymm7

        vpmullw 0x20(%rcx),  %ymm4,  %ymm4
        vpmullw 0x20(%rcx),  %ymm5,  %ymm5
        vpmulhw 0x00(%rcx),  %ymm4,  %ymm4
        vpmulhw 0x00(%rcx),  %ymm5,  %ymm5

        vpsubw       %ymm4,  %ymm6,  %ymm6
        vpsubw       %ymm5,  %ymm7,  %ymm7

        vpsubw       %ymm7,  %ymm6,  %ymm6

        vpaddw       %ymm6,  %ymm2,  %ymm4
        vpsubw       %ymm6,  %ymm2,  %ymm5

        vmovdqu      %ymm4, 0x00(%rdi)
        vmovdqu      %ymm5, 0x40(%rdi)

// ================

        vpmullw     %ymm15, %ymm10,  %ymm4
        vpmulhw     %ymm15, %ymm10,  %ymm6
        vpmullw     %ymm14, %ymm11,  %ymm5
        vpmulhw     %ymm14, %ymm11,  %ymm7

        vpmullw 0x20(%rcx),  %ymm4,  %ymm4
        vpmullw 0x20(%rcx),  %ymm5,  %ymm5
        vpmulhw 0x00(%rcx),  %ymm4,  %ymm4
        vpmulhw 0x00(%rcx),  %ymm5,  %ymm5

        vpsubw       %ymm4,  %ymm6,  %ymm6
        vpsubw       %ymm5,  %ymm7,  %ymm7

        vpaddw       %ymm7,  %ymm6,  %ymm6

        vpaddw       %ymm6,  %ymm3,  %ymm4
        vpsubw       %ymm6,  %ymm3,  %ymm5

        vmovdqu      %ymm4, 0x20(%rdi)
        vmovdqu      %ymm5, 0x60(%rdi)


// ================================

        ret

// x^4 + 1 =
// (x^2 + sqrt(2) x + 1) (x^2 - sqrt(2) x + 1)

// (a0, a1, a2, a3) =>
// (a0 - a2, a1 + a3, sqrt(2) a2, sqrt(2) a3) =>
// (a0 - a2, a1 + a3, sqrt(2) a3, sqrt(2) a2) =>
// (a0 - a2 + sqrt(2) a3, a1 + a3 - sqrt(2) a3, a0 - a2 - sqrt(2) a3 , a1 + a3 + sqrt(2) a3)

.global __asm_negacyclic_FFT4
.global ___asm_negacyclic_FFT4
__asm_negacyclic_FFT4:
___asm_negacyclic_FFT4:

        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x080, %rsp

        // a0, a1, a2, a3
        vmovdqu 0x00(%rsi), %ymm0
        vmovdqu 0x20(%rsi), %ymm1
        vmovdqu 0x40(%rsi), %ymm2
        vmovdqu 0x60(%rsi), %ymm3

        // b0, b1, b2, b3

        vmovdqu 0x00(%rdx), %ymm4
        vmovdqu 0x20(%rdx), %ymm5
        vmovdqu 0x40(%rdx), %ymm6
        vmovdqu 0x60(%rdx), %ymm7

// ================================

        vpsubw       %ymm2,  %ymm0,  %ymm0
        vpaddw       %ymm3,  %ymm1,  %ymm1

        vpsubw       %ymm6,  %ymm4,  %ymm4
        vpaddw       %ymm7,  %ymm5,  %ymm5

        montgomery_mul_mem 2, 2, 0x060(%rcx), 0x000(%rcx), 0x020(%rcx), 14, 15
        montgomery_mul_mem 3, 3, 0x060(%rcx), 0x000(%rcx), 0x020(%rcx), 14, 15
        montgomery_mul_mem 6, 6, 0x060(%rcx), 0x000(%rcx), 0x020(%rcx), 14, 15
        montgomery_mul_mem 7, 7, 0x060(%rcx), 0x000(%rcx), 0x020(%rcx), 14, 15

        vpaddw       %ymm7,  %ymm4,  %ymm8
        vpsubw       %ymm7,  %ymm4, %ymm10
        vpsubw       %ymm6,  %ymm5,  %ymm9
        vpaddw       %ymm6,  %ymm5, %ymm11

        vpaddw       %ymm3,  %ymm0,  %ymm4
        vpsubw       %ymm3,  %ymm0,  %ymm6
        vpsubw       %ymm2,  %ymm1,  %ymm5
        vpaddw       %ymm2,  %ymm1,  %ymm7

// ================

        montgomery_mul_constmem  1,  4,  9, 0x000(%rcx), 0x020(%rcx), 14, 15
        montgomery_mul_constmem  2,  5,  8, 0x000(%rcx), 0x020(%rcx), 14, 15
        vpaddw       %ymm2,  %ymm1,  %ymm1
        montgomery_mul_constmem  0,  4,  8, 0x000(%rcx), 0x020(%rcx), 14, 15
        montgomery_mul_constmem  2,  5,  9, 0x000(%rcx), 0x020(%rcx), 14, 15
        vpsubw       %ymm2,  %ymm0,  %ymm0
        vmovdqu 0x60(%rcx), %ymm3
        montgomery_mul_constmem  2,  2,  3, 0x000(%rcx), 0x020(%rcx), 14, 15
        vpsubw       %ymm2,  %ymm1,  %ymm1

        montgomery_mul_constmem  3,  6, 11, 0x000(%rcx), 0x020(%rcx), 14, 15
        montgomery_mul_constmem  4,  7, 10, 0x000(%rcx), 0x020(%rcx), 14, 15
        vpaddw       %ymm4,  %ymm3,  %ymm3
        montgomery_mul_constmem  2,  6, 10, 0x000(%rcx), 0x020(%rcx), 14, 15
        montgomery_mul_constmem  4,  7, 11, 0x000(%rcx), 0x020(%rcx), 14, 15
        vpsubw       %ymm4,  %ymm2,  %ymm2
        vmovdqu 0x60(%rcx), %ymm11
        montgomery_mul_constmem  4,  4, 11, 0x000(%rcx), 0x020(%rcx), 14, 15
        vpaddw       %ymm4,  %ymm3,  %ymm3

// ================================

        vpaddw       %ymm2,  %ymm0,  %ymm4
        vpsubw       %ymm2,  %ymm0,  %ymm6
        vpaddw       %ymm1,  %ymm3,  %ymm5
        vpsubw       %ymm1,  %ymm3,  %ymm7

        montgomery_mul_mem 6, 6, 0x080(%rcx), 0x000(%rcx), 0x020(%rcx), 14, 15
        montgomery_mul_mem 7, 7, 0x080(%rcx), 0x000(%rcx), 0x020(%rcx), 14, 15


        vpaddw       %ymm7,  %ymm4,  %ymm4
        vpsubw       %ymm6,  %ymm5,  %ymm5


        vmovdqu      %ymm4, 0x00(%rdi)
        vmovdqu      %ymm5, 0x20(%rdi)
        vmovdqu      %ymm7, 0x40(%rdi)
        vmovdqu      %ymm6, 0x60(%rdi)

        vzeroupper
        leave
        ret

.global __asm_cyclic_FFT8
.global ___asm_cyclic_FFT8
__asm_cyclic_FFT8:
___asm_cyclic_FFT8:

        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x100, %rsp

        vmovdqu 0x00(%rsi), %ymm0
        vmovdqu 0x20(%rsi), %ymm1
        vmovdqu 0x40(%rsi), %ymm2
        vmovdqu 0x60(%rsi), %ymm3
        vmovdqu 0x80(%rsi), %ymm4
        vmovdqu 0xa0(%rsi), %ymm5
        vmovdqu 0xc0(%rsi), %ymm6
        vmovdqu 0xe0(%rsi), %ymm7

        vpaddw       %ymm4,  %ymm0,  %ymm8
        vpsubw       %ymm4,  %ymm0, %ymm12
        vpaddw       %ymm5,  %ymm1,  %ymm9
        vpsubw       %ymm5,  %ymm1, %ymm13
        vpaddw       %ymm6,  %ymm2, %ymm10
        vpsubw       %ymm6,  %ymm2, %ymm14
        vpaddw       %ymm7,  %ymm3, %ymm11
        vpsubw       %ymm7,  %ymm3, %ymm15

        vpaddw      %ymm10,  %ymm8,  %ymm0
        vpsubw      %ymm10,  %ymm8,  %ymm2
        vpaddw      %ymm11,  %ymm9,  %ymm1
        vpsubw      %ymm11,  %ymm9,  %ymm3

        vmovdqu      %ymm0, 0x00(%rsp)
        vmovdqu      %ymm1, 0x20(%rsp)
        vmovdqu      %ymm2, 0x40(%rsp)
        vmovdqu      %ymm3, 0x60(%rsp)
        vmovdqu     %ymm12, 0x80(%rsp)
        vmovdqu     %ymm13, 0xa0(%rsp)
        vmovdqu     %ymm14, 0xc0(%rsp)
        vmovdqu     %ymm15, 0xe0(%rsp)

        vmovdqu 0x00(%rdx), %ymm0
        vmovdqu 0x20(%rdx), %ymm1
        vmovdqu 0x40(%rdx), %ymm2
        vmovdqu 0x60(%rdx), %ymm3
        vmovdqu 0x80(%rdx), %ymm4
        vmovdqu 0xa0(%rdx), %ymm5
        vmovdqu 0xc0(%rdx), %ymm6
        vmovdqu 0xe0(%rdx), %ymm7

        vpaddw       %ymm4,  %ymm0,  %ymm8
        vpsubw       %ymm4,  %ymm0, %ymm12
        vpaddw       %ymm5,  %ymm1,  %ymm9
        vpsubw       %ymm5,  %ymm1, %ymm13
        vpaddw       %ymm6,  %ymm2, %ymm10
        vpsubw       %ymm6,  %ymm2, %ymm14
        vpaddw       %ymm7,  %ymm3, %ymm11
        vpsubw       %ymm7,  %ymm3, %ymm15

        vpaddw      %ymm10,  %ymm8,  %ymm0
        vpsubw      %ymm10,  %ymm8,  %ymm2
        vpaddw      %ymm11,  %ymm9,  %ymm1
        vpsubw      %ymm11,  %ymm9,  %ymm3

        vmovdqu 0x00(%rsp),  %ymm8
        vmovdqu 0x20(%rsp),  %ymm9
        vmovdqu 0x40(%rsp), %ymm10
        vmovdqu 0x60(%rsp), %ymm11

        vmovdqu     %ymm12, 0x00(%rsp)
        vmovdqu     %ymm13, 0x20(%rsp)
        vmovdqu     %ymm14, 0x40(%rsp)
        vmovdqu     %ymm15, 0x60(%rsp)

        vmovdqu 0x00(%rcx),  %ymm12
        vmovdqu 0x20(%rcx),  %ymm13

        __asm_cyclic_schoolbook2_core %rdi, 0x00, \
                0, 1, 8, 9, \
                12, 13, \
                4, 5, 6, 7, 14, 15

        __asm_negacyclic_schoolbook2_core %rdi, 0x40, \
                2, 3, 10, 11, \
                12, 13, \
                4, 5, 6, 7, 14, 15

        vmovdqu 0x80(%rsp),  %ymm0
        vmovdqu 0xa0(%rsp),  %ymm1
        vmovdqu 0xc0(%rsp),  %ymm2
        vmovdqu 0xe0(%rsp),  %ymm3
        vmovdqu 0x00(%rsp),  %ymm8
        vmovdqu 0x20(%rsp),  %ymm9
        vmovdqu 0x40(%rsp), %ymm10
        vmovdqu 0x60(%rsp), %ymm11

        __asm_negacyclic_schoolbook4_core %rdi, 0x80, \
                0, 1, 2, 3, \
                8, 9, 10, 11, \
                12, 13, 0x20(%rcx), 0x40(%rcx), \
                4, 5, 6, 7, 14, 15

        vmovdqu 0x00(%rdi),  %ymm0
        vmovdqu 0x20(%rdi),  %ymm1
        vmovdqu 0x40(%rdi),  %ymm2
        vmovdqu 0x60(%rdi),  %ymm3
        vmovdqu 0x80(%rdi),  %ymm4
        vmovdqu 0xa0(%rdi),  %ymm5
        vmovdqu 0xc0(%rdi),  %ymm6
        vmovdqu 0xe0(%rdi),  %ymm7

        vpaddw       %ymm2,  %ymm0,  %ymm8
        vpsubw       %ymm2,  %ymm0, %ymm10
        vpaddw       %ymm3,  %ymm1,  %ymm9
        vpsubw       %ymm3,  %ymm1, %ymm11

        vmovdqu      %ymm8,  %ymm0
        vmovdqu      %ymm9,  %ymm1
        vmovdqu     %ymm10,  %ymm2
        vmovdqu     %ymm11,  %ymm3

        vpaddw       %ymm4,  %ymm4,  %ymm4
        vpaddw       %ymm5,  %ymm5,  %ymm5
        vpaddw       %ymm6,  %ymm6,  %ymm6
        vpaddw       %ymm7,  %ymm7,  %ymm7

        vpaddw       %ymm4,  %ymm0,  %ymm8
        vpsubw       %ymm4,  %ymm0, %ymm12
        vpaddw       %ymm5,  %ymm1,  %ymm9
        vpsubw       %ymm5,  %ymm1, %ymm13
        vpaddw       %ymm6,  %ymm2, %ymm10
        vpsubw       %ymm6,  %ymm2, %ymm14
        vpaddw       %ymm7,  %ymm3, %ymm11
        vpsubw       %ymm7,  %ymm3, %ymm15

        vmovdqu      %ymm8, 0x00(%rdi)
        vmovdqu      %ymm9, 0x20(%rdi)
        vmovdqu     %ymm10, 0x40(%rdi)
        vmovdqu     %ymm11, 0x60(%rdi)
        vmovdqu     %ymm12, 0x80(%rdi)
        vmovdqu     %ymm13, 0xa0(%rdi)
        vmovdqu     %ymm14, 0xc0(%rdi)
        vmovdqu     %ymm15, 0xe0(%rdi)

        vzeroupper
        leave
        ret

.global __asm_karatsuba8
.global ___asm_karatsuba8
__asm_karatsuba8:
___asm_karatsuba8:


        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x300, %rsp

        __asm_karatsuba8_core \
                %rdi, %rsi, %rdx, %rsp, 0x00, 0x00, 0x00, 0x00, \
                0x00(%rcx), 0x20(%rcx), 0x40(%rcx)

        vzeroupper
        leave
        ret

.global __asm_negacyclic_karatsuba8
.global ___asm_negacyclic_karatsuba8
__asm_negacyclic_karatsuba8:
___asm_negacyclic_karatsuba8:

        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x300, %rsp

        __asm_negacyclic_karatsuba8_core \
                %rdi, %rsi, %rdx, %rsp, 0x00, 0x00, 0x00, 0x00, \
                0x00(%rcx), 0x20(%rcx), 0x40(%rcx)

        vzeroupper
        leave
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

.global __asm_cyclic_FFT16_precompute
.global ___asm_cyclic_FFT16_precompute
__asm_cyclic_FFT16_precompute:
___asm_cyclic_FFT16_precompute:

        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x600, %rsp

        vmovdqu 0x000(%rsi),  %ymm0
        vmovdqu 0x020(%rsi),  %ymm1
        vmovdqu 0x040(%rsi),  %ymm2
        vmovdqu 0x060(%rsi),  %ymm3
        vmovdqu 0x080(%rsi),  %ymm4
        vmovdqu 0x0a0(%rsi),  %ymm5
        vmovdqu 0x0c0(%rsi),  %ymm6
        vmovdqu 0x0e0(%rsi),  %ymm7

        vpsubw  0x100(%rsi),  %ymm0,  %ymm8
        vpaddw  0x100(%rsi),  %ymm0,  %ymm0
        vpsubw  0x120(%rsi),  %ymm1,  %ymm9
        vpaddw  0x120(%rsi),  %ymm1,  %ymm1
        vpsubw  0x140(%rsi),  %ymm2, %ymm10
        vpaddw  0x140(%rsi),  %ymm2,  %ymm2
        vpsubw  0x160(%rsi),  %ymm3, %ymm11
        vpaddw  0x160(%rsi),  %ymm3,  %ymm3
        vpsubw  0x180(%rsi),  %ymm4, %ymm12
        vpaddw  0x180(%rsi),  %ymm4,  %ymm4
        vpsubw  0x1a0(%rsi),  %ymm5, %ymm13
        vpaddw  0x1a0(%rsi),  %ymm5,  %ymm5
        vpsubw  0x1c0(%rsi),  %ymm6, %ymm14
        vpaddw  0x1c0(%rsi),  %ymm6,  %ymm6
        vpsubw  0x1e0(%rsi),  %ymm7, %ymm15
        vpaddw  0x1e0(%rsi),  %ymm7,  %ymm7

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

        vmovdqu       %ymm0, 0x000(%rdi)
        vmovdqu       %ymm1, 0x020(%rdi)
        vmovdqu       %ymm2, 0x040(%rdi)
        vmovdqu       %ymm3, 0x060(%rdi)
        vmovdqu       %ymm4, 0x080(%rdi)
        vmovdqu       %ymm5, 0x0a0(%rdi)
        vmovdqu       %ymm6, 0x0c0(%rdi)
        vmovdqu       %ymm7, 0x0e0(%rdi)
        vmovdqu       %ymm8, 0x100(%rdi)
        vmovdqu       %ymm9, 0x120(%rdi)
        vmovdqu      %ymm10, 0x140(%rdi)
        vmovdqu      %ymm11, 0x160(%rdi)
        vmovdqu      %ymm12, 0x180(%rdi)
        vmovdqu      %ymm13, 0x1a0(%rdi)
        vmovdqu      %ymm14, 0x1c0(%rdi)
        vmovdqu      %ymm15, 0x1e0(%rdi)



        vzeroupper
        leave
        ret

.if 0

src1mid

0x000
0x020
0x040
0x060
0x080
0x0a0
0x0c0
0x0e0

src2mid

0x100
0x120
0x140
0x160
0x180
0x1a0
0x1c0
0x1e0

desmid

0x200
0x220
0x240
0x260
0x280
0x2a0
0x2c0
0x2e0
0x300
0x320
0x340
0x360
0x380
0x3a0
0x3c0
0x3e0

deshi

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

.global __asm_weighted_karatsuba16
.global ___asm_weighted_karatsuba16
__asm_weighted_karatsuba16:
___asm_weighted_karatsuba16:

        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x900, %rsp

        vmovdqu 0x000(%rsi), %ymm0
        vmovdqu 0x020(%rsi), %ymm1
        vmovdqu 0x040(%rsi), %ymm2
        vmovdqu 0x060(%rsi), %ymm3
        vmovdqu 0x080(%rsi), %ymm4
        vmovdqu 0x0a0(%rsi), %ymm5
        vmovdqu 0x0c0(%rsi), %ymm6
        vmovdqu 0x0e0(%rsi), %ymm7

        vpaddw  0x100(%rsi), %ymm0, %ymm0
        vpaddw  0x120(%rsi), %ymm1, %ymm1
        vpaddw  0x140(%rsi), %ymm2, %ymm2
        vpaddw  0x160(%rsi), %ymm3, %ymm3
        vpaddw  0x180(%rsi), %ymm4, %ymm4
        vpaddw  0x1a0(%rsi), %ymm5, %ymm5
        vpaddw  0x1c0(%rsi), %ymm6, %ymm6
        vpaddw  0x1e0(%rsi), %ymm7, %ymm7

        vmovdqu       %ymm0, 0x000(%rsp)
        vmovdqu       %ymm1, 0x020(%rsp)
        vmovdqu       %ymm2, 0x040(%rsp)
        vmovdqu       %ymm3, 0x060(%rsp)
        vmovdqu       %ymm4, 0x080(%rsp)
        vmovdqu       %ymm5, 0x0a0(%rsp)
        vmovdqu       %ymm6, 0x0c0(%rsp)
        vmovdqu       %ymm7, 0x0e0(%rsp)

        vmovdqu 0x000(%rdx), %ymm0
        vmovdqu 0x020(%rdx), %ymm1
        vmovdqu 0x040(%rdx), %ymm2
        vmovdqu 0x060(%rdx), %ymm3
        vmovdqu 0x080(%rdx), %ymm4
        vmovdqu 0x0a0(%rdx), %ymm5
        vmovdqu 0x0c0(%rdx), %ymm6
        vmovdqu 0x0e0(%rdx), %ymm7

        vpaddw  0x100(%rdx), %ymm0, %ymm0
        vpaddw  0x120(%rdx), %ymm1, %ymm1
        vpaddw  0x140(%rdx), %ymm2, %ymm2
        vpaddw  0x160(%rdx), %ymm3, %ymm3
        vpaddw  0x180(%rdx), %ymm4, %ymm4
        vpaddw  0x1a0(%rdx), %ymm5, %ymm5
        vpaddw  0x1c0(%rdx), %ymm6, %ymm6
        vpaddw  0x1e0(%rdx), %ymm7, %ymm7

        vmovdqu       %ymm0, 0x100(%rsp)
        vmovdqu       %ymm1, 0x120(%rsp)
        vmovdqu       %ymm2, 0x140(%rsp)
        vmovdqu       %ymm3, 0x160(%rsp)
        vmovdqu       %ymm4, 0x180(%rsp)
        vmovdqu       %ymm5, 0x1a0(%rsp)
        vmovdqu       %ymm6, 0x1c0(%rsp)
        vmovdqu       %ymm7, 0x1e0(%rsp)

        __asm_karatsuba8_core \
                %rsp, %rsp, %rsp, %rsp, 0x200, 0x000, 0x100, 0x600, \
                0x00(%rcx), 0x20(%rcx), 0x40(%rcx)

        __asm_karatsuba8_core \
                %rdi, %rsi, %rdx, %rsp, 0x000, 0x000, 0x000, 0x600, \
                0x00(%rcx), 0x20(%rcx), 0x40(%rcx)

        __asm_karatsuba8_core \
                %rsp, %rsi, %rdx, %rsp, 0x000, 0x100, 0x100, 0x600, \
                0x00(%rcx), 0x20(%rcx), 0x40(%rcx)


        vmovdqu 0x200(%rsp),  %ymm0
        vmovdqu 0x220(%rsp),  %ymm1
        vmovdqu 0x240(%rsp),  %ymm2
        vmovdqu 0x260(%rsp),  %ymm3
        vmovdqu 0x280(%rsp),  %ymm4
        vmovdqu 0x2a0(%rsp),  %ymm5
        vmovdqu 0x2c0(%rsp),  %ymm6
        vmovdqu 0x2e0(%rsp),  %ymm7
        vmovdqu 0x300(%rsp),  %ymm8
        vmovdqu 0x320(%rsp),  %ymm9
        vmovdqu 0x340(%rsp), %ymm10
        vmovdqu 0x360(%rsp), %ymm11
        vmovdqu 0x380(%rsp), %ymm12
        vmovdqu 0x3a0(%rsp), %ymm13
        vmovdqu 0x3c0(%rsp), %ymm14

        vpsubw  0x000(%rdi),  %ymm0,  %ymm0
        vpsubw  0x020(%rdi),  %ymm1,  %ymm1
        vpsubw  0x040(%rdi),  %ymm2,  %ymm2
        vpsubw  0x060(%rdi),  %ymm3,  %ymm3
        vpsubw  0x080(%rdi),  %ymm4,  %ymm4
        vpsubw  0x0a0(%rdi),  %ymm5,  %ymm5
        vpsubw  0x0c0(%rdi),  %ymm6,  %ymm6
        vpsubw  0x0e0(%rdi),  %ymm7,  %ymm7
        vpsubw  0x100(%rdi),  %ymm8,  %ymm8
        vpsubw  0x120(%rdi),  %ymm9,  %ymm9
        vpsubw  0x140(%rdi), %ymm10, %ymm10
        vpsubw  0x160(%rdi), %ymm11, %ymm11
        vpsubw  0x180(%rdi), %ymm12, %ymm12
        vpsubw  0x1a0(%rdi), %ymm13, %ymm13
        vpsubw  0x1c0(%rdi), %ymm14, %ymm14

        vpsubw  0x000(%rsp),  %ymm0,  %ymm0
        vpsubw  0x020(%rsp),  %ymm1,  %ymm1
        vpsubw  0x040(%rsp),  %ymm2,  %ymm2
        vpsubw  0x060(%rsp),  %ymm3,  %ymm3
        vpsubw  0x080(%rsp),  %ymm4,  %ymm4
        vpsubw  0x0a0(%rsp),  %ymm5,  %ymm5
        vpsubw  0x0c0(%rsp),  %ymm6,  %ymm6
        vpsubw  0x0e0(%rsp),  %ymm7,  %ymm7
        vpsubw  0x100(%rsp),  %ymm8,  %ymm8
        vpsubw  0x120(%rsp),  %ymm9,  %ymm9
        vpsubw  0x140(%rsp), %ymm10, %ymm10
        vpsubw  0x160(%rsp), %ymm11, %ymm11
        vpsubw  0x180(%rsp), %ymm12, %ymm12
        vpsubw  0x1a0(%rsp), %ymm13, %ymm13
        vpsubw  0x1c0(%rsp), %ymm14, %ymm14

        vmovdqu       %ymm7, 0x1e0(%rdi)

        vpaddw  0x000(%rsp),  %ymm8,  %ymm8
        vpaddw  0x020(%rsp),  %ymm9,  %ymm9
        vpaddw  0x040(%rsp), %ymm10, %ymm10
        vpaddw  0x060(%rsp), %ymm11, %ymm11
        vpaddw  0x080(%rsp), %ymm12, %ymm12
        vpaddw  0x0a0(%rsp), %ymm13, %ymm13
        vpaddw  0x0c0(%rsp), %ymm14, %ymm14

        vmovdqu       %ymm8, 0x000(%rsp)
        vmovdqu       %ymm9, 0x020(%rsp)
        vmovdqu      %ymm10, 0x040(%rsp)
        vmovdqu      %ymm11, 0x060(%rsp)
        vmovdqu      %ymm12, 0x080(%rsp)
        vmovdqu      %ymm13, 0x0a0(%rsp)
        vmovdqu      %ymm14, 0x0c0(%rsp)

// reduction: 0 ~ 6

        vmovdqu 0x000( %r8), %ymm15
        vpmullw 0x020(%rcx), %ymm15,  %ymm7

        vpmullw 0x0e0(%rsp),  %ymm7,  %ymm8
        vpmulhw 0x0e0(%rsp), %ymm15, %ymm12
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpsubw        %ymm8, %ymm12, %ymm12
        vpaddw  0x0e0(%rdi), %ymm12, %ymm12
        vmovdqu      %ymm12, 0x0e0(%rdi)

        vpmullw 0x000(%rsp),  %ymm7,  %ymm8
        vpmullw 0x020(%rsp),  %ymm7,  %ymm9
        vpmullw 0x040(%rsp),  %ymm7, %ymm10
        vpmulhw 0x000(%rsp), %ymm15, %ymm12
        vpmulhw 0x020(%rsp), %ymm15, %ymm13
        vpmulhw 0x040(%rsp), %ymm15, %ymm14
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpmulhw 0x000(%rcx),  %ymm9,  %ymm9
        vpmulhw 0x000(%rcx), %ymm10, %ymm10
        vpsubw        %ymm8, %ymm12, %ymm12
        vpsubw        %ymm9, %ymm13, %ymm13
        vpsubw       %ymm10, %ymm14, %ymm14
        vmovdqu 0x000(%rdi),  %ymm8
        vmovdqu 0x020(%rdi),  %ymm9
        vmovdqu 0x040(%rdi), %ymm10
        vpaddw       %ymm12,  %ymm8,  %ymm8
        vpaddw       %ymm13,  %ymm9,  %ymm9
        vpaddw       %ymm14, %ymm10, %ymm10
        vmovdqu       %ymm8, 0x000(%rdi)
        vmovdqu       %ymm9, 0x020(%rdi)
        vmovdqu      %ymm10, 0x040(%rdi)

        vpmullw 0x060(%rsp),  %ymm7,  %ymm8
        vpmullw 0x080(%rsp),  %ymm7,  %ymm9
        vpmullw 0x0a0(%rsp),  %ymm7, %ymm10
        vpmullw 0x0c0(%rsp),  %ymm7, %ymm11
        vpmulhw 0x060(%rsp), %ymm15, %ymm12
        vpmulhw 0x080(%rsp), %ymm15, %ymm13
        vpmulhw 0x0a0(%rsp), %ymm15, %ymm14
        vpmulhw 0x0c0(%rsp), %ymm15, %ymm15
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpmulhw 0x000(%rcx),  %ymm9,  %ymm9
        vpmulhw 0x000(%rcx), %ymm10, %ymm10
        vpmulhw 0x000(%rcx), %ymm11, %ymm11
        vpsubw        %ymm8, %ymm12, %ymm12
        vpsubw        %ymm9, %ymm13, %ymm13
        vpsubw       %ymm10, %ymm14, %ymm14
        vpsubw       %ymm11, %ymm15, %ymm15
        vmovdqu 0x060(%rdi),  %ymm8
        vmovdqu 0x080(%rdi),  %ymm9
        vmovdqu 0x0a0(%rdi), %ymm10
        vmovdqu 0x0c0(%rdi), %ymm11
        vpaddw       %ymm12,  %ymm8,  %ymm8
        vpaddw       %ymm13,  %ymm9,  %ymm9
        vpaddw       %ymm14, %ymm10, %ymm10
        vpaddw       %ymm15, %ymm11, %ymm11
        vmovdqu       %ymm8, 0x060(%rdi)
        vmovdqu       %ymm9, 0x080(%rdi)
        vmovdqu      %ymm10, 0x0a0(%rdi)
        vmovdqu      %ymm11, 0x0c0(%rdi)

        vpaddw  0x100(%rdi),  %ymm0,  %ymm0
        vpaddw  0x120(%rdi),  %ymm1,  %ymm1
        vpaddw  0x140(%rdi),  %ymm2,  %ymm2
        vpaddw  0x160(%rdi),  %ymm3,  %ymm3
        vpaddw  0x180(%rdi),  %ymm4,  %ymm4
        vpaddw  0x1a0(%rdi),  %ymm5,  %ymm5
        vpaddw  0x1c0(%rdi),  %ymm6,  %ymm6

// reduction 8 ~ 14

        vmovdqu 0x000( %r8), %ymm15
        vpmullw 0x020(%rcx), %ymm15,  %ymm7

        vpmullw 0x100(%rsp),  %ymm7,  %ymm8
        vpmullw 0x120(%rsp),  %ymm7,  %ymm9
        vpmullw 0x140(%rsp),  %ymm7, %ymm10
        vpmulhw 0x100(%rsp), %ymm15, %ymm12
        vpmulhw 0x120(%rsp), %ymm15, %ymm13
        vpmulhw 0x140(%rsp), %ymm15, %ymm14
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpmulhw 0x000(%rcx),  %ymm9,  %ymm9
        vpmulhw 0x000(%rcx), %ymm10, %ymm10
        vpsubw        %ymm8, %ymm12, %ymm12
        vpsubw        %ymm9, %ymm13, %ymm13
        vpsubw       %ymm10, %ymm14, %ymm14
        vpaddw        %ymm0, %ymm12,  %ymm0
        vpaddw        %ymm1, %ymm13,  %ymm1
        vpaddw        %ymm2, %ymm14,  %ymm2

        vpmullw 0x160(%rsp),  %ymm7,  %ymm8
        vpmullw 0x180(%rsp),  %ymm7,  %ymm9
        vpmullw 0x1a0(%rsp),  %ymm7, %ymm10
        vpmullw 0x1c0(%rsp),  %ymm7, %ymm11
        vpmulhw 0x160(%rsp), %ymm15, %ymm12
        vpmulhw 0x180(%rsp), %ymm15, %ymm13
        vpmulhw 0x1a0(%rsp), %ymm15, %ymm14
        vpmulhw 0x1c0(%rsp), %ymm15, %ymm15
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpmulhw 0x000(%rcx),  %ymm9,  %ymm9
        vpmulhw 0x000(%rcx), %ymm10, %ymm10
        vpmulhw 0x000(%rcx), %ymm11, %ymm11
        vpsubw        %ymm8, %ymm12, %ymm12
        vpsubw        %ymm9, %ymm13, %ymm13
        vpsubw       %ymm10, %ymm14, %ymm14
        vpsubw       %ymm11, %ymm15, %ymm15
        vpaddw        %ymm3, %ymm12,  %ymm3
        vpaddw        %ymm4, %ymm13,  %ymm4
        vpaddw        %ymm5, %ymm14,  %ymm5
        vpaddw        %ymm6, %ymm15,  %ymm6

        vmovdqu       %ymm0, 0x100(%rdi)
        vmovdqu       %ymm1, 0x120(%rdi)
        vmovdqu       %ymm2, 0x140(%rdi)
        vmovdqu       %ymm3, 0x160(%rdi)
        vmovdqu       %ymm4, 0x180(%rdi)
        vmovdqu       %ymm5, 0x1a0(%rdi)
        vmovdqu       %ymm6, 0x1c0(%rdi)

        vzeroupper
        leave
        ret


.global __asm_weighted_double_karatsuba16
.global ___asm_weighted_double_karatsuba16
__asm_weighted_double_karatsuba16:
___asm_weighted_double_karatsuba16:

        push   %rbp
        mov    %rsp, %rbp
        and    $0xffffffffffffffe0, %rsp
        sub    $0x900, %rsp

        vmovdqu 0x000(%rsi), %ymm0
        vmovdqu 0x020(%rsi), %ymm1
        vmovdqu 0x040(%rsi), %ymm2
        vmovdqu 0x060(%rsi), %ymm3
        vmovdqu 0x080(%rsi), %ymm4
        vmovdqu 0x0a0(%rsi), %ymm5
        vmovdqu 0x0c0(%rsi), %ymm6
        vmovdqu 0x0e0(%rsi), %ymm7

        vpaddw  0x100(%rsi), %ymm0, %ymm0
        vpaddw  0x120(%rsi), %ymm1, %ymm1
        vpaddw  0x140(%rsi), %ymm2, %ymm2
        vpaddw  0x160(%rsi), %ymm3, %ymm3
        vpaddw  0x180(%rsi), %ymm4, %ymm4
        vpaddw  0x1a0(%rsi), %ymm5, %ymm5
        vpaddw  0x1c0(%rsi), %ymm6, %ymm6
        vpaddw  0x1e0(%rsi), %ymm7, %ymm7

        vmovdqu       %ymm0, 0x000(%rsp)
        vmovdqu       %ymm1, 0x020(%rsp)
        vmovdqu       %ymm2, 0x040(%rsp)
        vmovdqu       %ymm3, 0x060(%rsp)
        vmovdqu       %ymm4, 0x080(%rsp)
        vmovdqu       %ymm5, 0x0a0(%rsp)
        vmovdqu       %ymm6, 0x0c0(%rsp)
        vmovdqu       %ymm7, 0x0e0(%rsp)

        vmovdqu 0x000(%rdx), %ymm0
        vmovdqu 0x020(%rdx), %ymm1
        vmovdqu 0x040(%rdx), %ymm2
        vmovdqu 0x060(%rdx), %ymm3
        vmovdqu 0x080(%rdx), %ymm4
        vmovdqu 0x0a0(%rdx), %ymm5
        vmovdqu 0x0c0(%rdx), %ymm6
        vmovdqu 0x0e0(%rdx), %ymm7

        vpaddw  0x100(%rdx), %ymm0, %ymm0
        vpaddw  0x120(%rdx), %ymm1, %ymm1
        vpaddw  0x140(%rdx), %ymm2, %ymm2
        vpaddw  0x160(%rdx), %ymm3, %ymm3
        vpaddw  0x180(%rdx), %ymm4, %ymm4
        vpaddw  0x1a0(%rdx), %ymm5, %ymm5
        vpaddw  0x1c0(%rdx), %ymm6, %ymm6
        vpaddw  0x1e0(%rdx), %ymm7, %ymm7

        vmovdqu       %ymm0, 0x100(%rsp)
        vmovdqu       %ymm1, 0x120(%rsp)
        vmovdqu       %ymm2, 0x140(%rsp)
        vmovdqu       %ymm3, 0x160(%rsp)
        vmovdqu       %ymm4, 0x180(%rsp)
        vmovdqu       %ymm5, 0x1a0(%rsp)
        vmovdqu       %ymm6, 0x1c0(%rsp)
        vmovdqu       %ymm7, 0x1e0(%rsp)

        __asm_karatsuba8_core \
                %rsp, %rsp, %rsp, %rsp, 0x200, 0x000, 0x100, 0x600, \
                0x00(%rcx), 0x20(%rcx), 0x40(%rcx)

        __asm_karatsuba8_core \
                %rdi, %rsi, %rdx, %rsp, 0x000, 0x000, 0x000, 0x600, \
                0x00(%rcx), 0x20(%rcx), 0x40(%rcx)

        __asm_karatsuba8_core \
                %rsp, %rsi, %rdx, %rsp, 0x000, 0x100, 0x100, 0x600, \
                0x00(%rcx), 0x20(%rcx), 0x40(%rcx)


        vmovdqu 0x200(%rsp),  %ymm0
        vmovdqu 0x220(%rsp),  %ymm1
        vmovdqu 0x240(%rsp),  %ymm2
        vmovdqu 0x260(%rsp),  %ymm3
        vmovdqu 0x280(%rsp),  %ymm4
        vmovdqu 0x2a0(%rsp),  %ymm5
        vmovdqu 0x2c0(%rsp),  %ymm6
        vmovdqu 0x2e0(%rsp),  %ymm7
        vmovdqu 0x300(%rsp),  %ymm8
        vmovdqu 0x320(%rsp),  %ymm9
        vmovdqu 0x340(%rsp), %ymm10
        vmovdqu 0x360(%rsp), %ymm11
        vmovdqu 0x380(%rsp), %ymm12
        vmovdqu 0x3a0(%rsp), %ymm13
        vmovdqu 0x3c0(%rsp), %ymm14

        vpsubw  0x000(%rdi),  %ymm0,  %ymm0
        vpsubw  0x020(%rdi),  %ymm1,  %ymm1
        vpsubw  0x040(%rdi),  %ymm2,  %ymm2
        vpsubw  0x060(%rdi),  %ymm3,  %ymm3
        vpsubw  0x080(%rdi),  %ymm4,  %ymm4
        vpsubw  0x0a0(%rdi),  %ymm5,  %ymm5
        vpsubw  0x0c0(%rdi),  %ymm6,  %ymm6
        vpsubw  0x0e0(%rdi),  %ymm7,  %ymm7
        vpsubw  0x100(%rdi),  %ymm8,  %ymm8
        vpsubw  0x120(%rdi),  %ymm9,  %ymm9
        vpsubw  0x140(%rdi), %ymm10, %ymm10
        vpsubw  0x160(%rdi), %ymm11, %ymm11
        vpsubw  0x180(%rdi), %ymm12, %ymm12
        vpsubw  0x1a0(%rdi), %ymm13, %ymm13
        vpsubw  0x1c0(%rdi), %ymm14, %ymm14

        vpsubw  0x000(%rsp),  %ymm0,  %ymm0
        vpsubw  0x020(%rsp),  %ymm1,  %ymm1
        vpsubw  0x040(%rsp),  %ymm2,  %ymm2
        vpsubw  0x060(%rsp),  %ymm3,  %ymm3
        vpsubw  0x080(%rsp),  %ymm4,  %ymm4
        vpsubw  0x0a0(%rsp),  %ymm5,  %ymm5
        vpsubw  0x0c0(%rsp),  %ymm6,  %ymm6
        vpsubw  0x0e0(%rsp),  %ymm7,  %ymm7
        vpsubw  0x100(%rsp),  %ymm8,  %ymm8
        vpsubw  0x120(%rsp),  %ymm9,  %ymm9
        vpsubw  0x140(%rsp), %ymm10, %ymm10
        vpsubw  0x160(%rsp), %ymm11, %ymm11
        vpsubw  0x180(%rsp), %ymm12, %ymm12
        vpsubw  0x1a0(%rsp), %ymm13, %ymm13
        vpsubw  0x1c0(%rsp), %ymm14, %ymm14

        vpaddw        %ymm7,  %ymm7,  %ymm7
        vmovdqu       %ymm7, 0x1e0(%rdi)

        vpaddw  0x000(%rsp),  %ymm8,  %ymm8
        vpaddw  0x020(%rsp),  %ymm9,  %ymm9
        vpaddw  0x040(%rsp), %ymm10, %ymm10
        vpaddw  0x060(%rsp), %ymm11, %ymm11
        vpaddw  0x080(%rsp), %ymm12, %ymm12
        vpaddw  0x0a0(%rsp), %ymm13, %ymm13
        vpaddw  0x0c0(%rsp), %ymm14, %ymm14

        vmovdqu       %ymm8, 0x000(%rsp)
        vmovdqu       %ymm9, 0x020(%rsp)
        vmovdqu      %ymm10, 0x040(%rsp)
        vmovdqu      %ymm11, 0x060(%rsp)
        vmovdqu      %ymm12, 0x080(%rsp)
        vmovdqu      %ymm13, 0x0a0(%rsp)
        vmovdqu      %ymm14, 0x0c0(%rsp)

// reduction: 0 ~ 6

        vmovdqu 0x000( %r8), %ymm15
        vpmullw 0x020(%rcx), %ymm15,  %ymm7

        vpmullw 0x0e0(%rsp),  %ymm7,  %ymm8
        vpmulhw 0x0e0(%rsp), %ymm15, %ymm12
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpsubw        %ymm8, %ymm12, %ymm12
        vpaddw  0x0e0(%rdi), %ymm12, %ymm12
        vpaddw       %ymm12, %ymm12, %ymm12
        barrett_reduce_constmem 12, 12, 0x000(%rcx), 0x040(%rcx), 8
        vmovdqu      %ymm12, 0x0e0(%rdi)

        vpmullw 0x000(%rsp),  %ymm7,  %ymm8
        vpmullw 0x020(%rsp),  %ymm7,  %ymm9
        vpmullw 0x040(%rsp),  %ymm7, %ymm10
        vpmulhw 0x000(%rsp), %ymm15, %ymm12
        vpmulhw 0x020(%rsp), %ymm15, %ymm13
        vpmulhw 0x040(%rsp), %ymm15, %ymm14
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpmulhw 0x000(%rcx),  %ymm9,  %ymm9
        vpmulhw 0x000(%rcx), %ymm10, %ymm10
        vpsubw        %ymm8, %ymm12, %ymm12
        vpsubw        %ymm9, %ymm13, %ymm13
        vpsubw       %ymm10, %ymm14, %ymm14
        vmovdqu 0x000(%rdi),  %ymm8
        vmovdqu 0x020(%rdi),  %ymm9
        vmovdqu 0x040(%rdi), %ymm10
        vpaddw       %ymm12,  %ymm8,  %ymm8
        vpaddw       %ymm13,  %ymm9,  %ymm9
        vpaddw       %ymm14, %ymm10, %ymm10
        vpaddw        %ymm8,  %ymm8,  %ymm8
        vpaddw        %ymm9,  %ymm9,  %ymm9
        vpaddw       %ymm10, %ymm10, %ymm10
        barrett_reduce_constmemx3 8, 9, 10, 8, 9, 10, 0x000(%rcx), 0x040(%rcx), 12, 13, 14
        vmovdqu       %ymm8, 0x000(%rdi)
        vmovdqu       %ymm9, 0x020(%rdi)
        vmovdqu      %ymm10, 0x040(%rdi)

        vpmullw 0x060(%rsp),  %ymm7,  %ymm8
        vpmullw 0x080(%rsp),  %ymm7,  %ymm9
        vpmullw 0x0a0(%rsp),  %ymm7, %ymm10
        vpmullw 0x0c0(%rsp),  %ymm7, %ymm11
        vpmulhw 0x060(%rsp), %ymm15, %ymm12
        vpmulhw 0x080(%rsp), %ymm15, %ymm13
        vpmulhw 0x0a0(%rsp), %ymm15, %ymm14
        vpmulhw 0x0c0(%rsp), %ymm15, %ymm15
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpmulhw 0x000(%rcx),  %ymm9,  %ymm9
        vpmulhw 0x000(%rcx), %ymm10, %ymm10
        vpmulhw 0x000(%rcx), %ymm11, %ymm11
        vpsubw        %ymm8, %ymm12, %ymm12
        vpsubw        %ymm9, %ymm13, %ymm13
        vpsubw       %ymm10, %ymm14, %ymm14
        vpsubw       %ymm11, %ymm15, %ymm15
        vmovdqu 0x060(%rdi),  %ymm8
        vmovdqu 0x080(%rdi),  %ymm9
        vmovdqu 0x0a0(%rdi), %ymm10
        vmovdqu 0x0c0(%rdi), %ymm11
        vpaddw       %ymm12,  %ymm8,  %ymm8
        vpaddw       %ymm13,  %ymm9,  %ymm9
        vpaddw       %ymm14, %ymm10, %ymm10
        vpaddw       %ymm15, %ymm11, %ymm11
        vpaddw        %ymm8,  %ymm8,  %ymm8
        vpaddw        %ymm9,  %ymm9,  %ymm9
        vpaddw       %ymm10, %ymm10, %ymm10
        vpaddw       %ymm11, %ymm11, %ymm11
        barrett_reduce_constmemx4 8, 9, 10, 11, 8, 9, 10, 11, 0x000(%rcx), 0x040(%rcx), 12, 13, 14, 15
        vmovdqu       %ymm8, 0x060(%rdi)
        vmovdqu       %ymm9, 0x080(%rdi)
        vmovdqu      %ymm10, 0x0a0(%rdi)
        vmovdqu      %ymm11, 0x0c0(%rdi)

        vpaddw  0x100(%rdi),  %ymm0,  %ymm0
        vpaddw  0x120(%rdi),  %ymm1,  %ymm1
        vpaddw  0x140(%rdi),  %ymm2,  %ymm2
        vpaddw  0x160(%rdi),  %ymm3,  %ymm3
        vpaddw  0x180(%rdi),  %ymm4,  %ymm4
        vpaddw  0x1a0(%rdi),  %ymm5,  %ymm5
        vpaddw  0x1c0(%rdi),  %ymm6,  %ymm6

// reduction 8 ~ 14

        vmovdqu 0x000( %r8), %ymm15
        vpmullw 0x020(%rcx), %ymm15,  %ymm7

        vpmullw 0x100(%rsp),  %ymm7,  %ymm8
        vpmullw 0x120(%rsp),  %ymm7,  %ymm9
        vpmullw 0x140(%rsp),  %ymm7, %ymm10
        vpmulhw 0x100(%rsp), %ymm15, %ymm12
        vpmulhw 0x120(%rsp), %ymm15, %ymm13
        vpmulhw 0x140(%rsp), %ymm15, %ymm14
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpmulhw 0x000(%rcx),  %ymm9,  %ymm9
        vpmulhw 0x000(%rcx), %ymm10, %ymm10
        vpsubw        %ymm8, %ymm12, %ymm12
        vpsubw        %ymm9, %ymm13, %ymm13
        vpsubw       %ymm10, %ymm14, %ymm14
        vpaddw        %ymm0, %ymm12,  %ymm0
        vpaddw        %ymm1, %ymm13,  %ymm1
        vpaddw        %ymm2, %ymm14,  %ymm2

        vpmullw 0x160(%rsp),  %ymm7,  %ymm8
        vpmullw 0x180(%rsp),  %ymm7,  %ymm9
        vpmullw 0x1a0(%rsp),  %ymm7, %ymm10
        vpmullw 0x1c0(%rsp),  %ymm7, %ymm11
        vpmulhw 0x160(%rsp), %ymm15, %ymm12
        vpmulhw 0x180(%rsp), %ymm15, %ymm13
        vpmulhw 0x1a0(%rsp), %ymm15, %ymm14
        vpmulhw 0x1c0(%rsp), %ymm15, %ymm15
        vpmulhw 0x000(%rcx),  %ymm8,  %ymm8
        vpmulhw 0x000(%rcx),  %ymm9,  %ymm9
        vpmulhw 0x000(%rcx), %ymm10, %ymm10
        vpmulhw 0x000(%rcx), %ymm11, %ymm11
        vpsubw        %ymm8, %ymm12, %ymm12
        vpsubw        %ymm9, %ymm13, %ymm13
        vpsubw       %ymm10, %ymm14, %ymm14
        vpsubw       %ymm11, %ymm15, %ymm15
        vpaddw        %ymm3, %ymm12,  %ymm3
        vpaddw        %ymm4, %ymm13,  %ymm4
        vpaddw        %ymm5, %ymm14,  %ymm5
        vpaddw        %ymm6, %ymm15,  %ymm6

        vpaddw        %ymm0,  %ymm0,  %ymm0
        vpaddw        %ymm1,  %ymm1,  %ymm1
        vpaddw        %ymm2,  %ymm2,  %ymm2
        vpaddw        %ymm3,  %ymm3,  %ymm3
        vpaddw        %ymm4,  %ymm4,  %ymm4
        vpaddw        %ymm5,  %ymm5,  %ymm5
        vpaddw        %ymm6,  %ymm6,  %ymm6

        barrett_reduce_constmemx4 0, 1, 2, 3, 0, 1, 2, 3, 0x000(%rcx), 0x040(%rcx), 8, 9, 10, 11
        barrett_reduce_constmemx3 4, 5, 6, 4, 5, 6, 0x000(%rcx), 0x040(%rcx), 12, 13, 14

        vmovdqu       %ymm0, 0x100(%rdi)
        vmovdqu       %ymm1, 0x120(%rdi)
        vmovdqu       %ymm2, 0x140(%rdi)
        vmovdqu       %ymm3, 0x160(%rdi)
        vmovdqu       %ymm4, 0x180(%rdi)
        vmovdqu       %ymm5, 0x1a0(%rdi)
        vmovdqu       %ymm6, 0x1c0(%rdi)


        vzeroupper
        leave
        ret

















