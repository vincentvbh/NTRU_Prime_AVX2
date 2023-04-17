
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "NTT_params.h"
#include "__avx2.h"


#include "cpucycles.h"

#define ITERATIONS 100000

uint64_t cycles[ITERATIONS + 1];
uint64_t t0, t1;

static int cmp_uint64(const void *a, const void *b){
    return (*(uint64_t*)a) - (*(uint64_t*)b);
}

int main(void){


    int16_t poly1[ARRAY_N], poly2[ARRAY_N];
    int16_t res[ARRAY_N];
    int16_t poly2inv[ARRAY_N];

    int16_t twiddle51[51];

    printf("8 repetitions\n");

// ================================
// 2 x 2

    printf("======== 2x2 ========\n");

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_schoolbook2(res, poly1, poly2, const_buff);
        __asm_schoolbook2(res, poly1, poly2, const_buff);
        __asm_schoolbook2(res, poly1, poly2, const_buff);
        __asm_schoolbook2(res, poly1, poly2, const_buff);
        __asm_schoolbook2(res, poly1, poly2, const_buff);
        __asm_schoolbook2(res, poly1, poly2, const_buff);
        __asm_schoolbook2(res, poly1, poly2, const_buff);
        __asm_schoolbook2(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm schoolbook 2x2: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm schoolbook cyclic 2x2: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm schoolbook negacyclic 2x2: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

// ================================
// 4 x 4

    printf("======== 4x4 ========\n");

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_schoolbook4(res, poly1, poly2, const_buff);
        __asm_schoolbook4(res, poly1, poly2, const_buff);
        __asm_schoolbook4(res, poly1, poly2, const_buff);
        __asm_schoolbook4(res, poly1, poly2, const_buff);
        __asm_schoolbook4(res, poly1, poly2, const_buff);
        __asm_schoolbook4(res, poly1, poly2, const_buff);
        __asm_schoolbook4(res, poly1, poly2, const_buff);
        __asm_schoolbook4(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm schoolbook 4x4: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm schoolbook cyclic 4x4: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
        __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm schoolbook negacyclic 4x4: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm FFT cyclic 4x4: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
        __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm FFT negacyclic 4x4: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

// ================================
// 8 x 8

    printf("======== 8x8 ========\n");

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_karatsuba8(res, poly1, poly2, const_buff);
        __asm_karatsuba8(res, poly1, poly2, const_buff);
        __asm_karatsuba8(res, poly1, poly2, const_buff);
        __asm_karatsuba8(res, poly1, poly2, const_buff);
        __asm_karatsuba8(res, poly1, poly2, const_buff);
        __asm_karatsuba8(res, poly1, poly2, const_buff);
        __asm_karatsuba8(res, poly1, poly2, const_buff);
        __asm_karatsuba8(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm Karatsuba 8x8: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
        __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm FFT cyclic 8x8: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
        __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
        __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
        __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
        __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
        __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
        __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
        __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm Karatsuba negacyclic 8x8: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

// ================================
// 16 x 16

    printf("======== 16x16 ========\n");

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
        __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
        __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
        __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
        __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
        __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
        __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
        __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm FFT cyclic 16x16: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_weighted_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm weighted 16x16: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_weighted_double_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_double_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_double_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_double_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_double_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_double_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_double_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        __asm_weighted_double_karatsuba16(res, poly1, poly2, const_buff, twiddle51);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm weighted doubling 16x16: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);
        weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);
        weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);
        weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);
        weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);
        weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);
        weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);
        weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("CT weighted 16x16: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);


// ================================
// butterflies

    printf("======== Butterflies ========\n");


    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_3x2(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm 3x2: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_3x2_pre(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_pre(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_pre(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_pre(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_pre(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_pre(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_pre(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_pre(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm 3x2 pre: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_3x2_post(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_post(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_post(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_post(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_post(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_post(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_post(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        __asm_3x2_post(poly1, poly1 + 17 * 48, OMEGA3_buff, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm 3x2 post: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);


    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_rader17(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm Rader-17: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        __asm_rader17_scalei(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17_scalei(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17_scalei(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17_scalei(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17_scalei(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17_scalei(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17_scalei(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        __asm_rader17_scalei(poly1, poly1, const_twiddle17, const_twiddle17inv, const_buff);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("asm Rader-17 scaled indices: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

// ================================

    printf("======== Transpose ========\n");


    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        transpose((int16_t (*)[16])poly1, (int16_t (*)[16])poly1);
        transpose((int16_t (*)[16])poly1, (int16_t (*)[16])poly1);
        transpose((int16_t (*)[16])poly1, (int16_t (*)[16])poly1);
        transpose((int16_t (*)[16])poly1, (int16_t (*)[16])poly1);
        transpose((int16_t (*)[16])poly1, (int16_t (*)[16])poly1);
        transpose((int16_t (*)[16])poly1, (int16_t (*)[16])poly1);
        transpose((int16_t (*)[16])poly1, (int16_t (*)[16])poly1);
        transpose((int16_t (*)[16])poly1, (int16_t (*)[16])poly1);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("transpose: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        transpose_partial3((int16_t (*)[16])poly2, poly1                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
        transpose_partial3((int16_t (*)[16])poly2, poly1                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
        transpose_partial3((int16_t (*)[16])poly2, poly1                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
        transpose_partial3((int16_t (*)[16])poly2, poly1                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
        transpose_partial3((int16_t (*)[16])poly2, poly1                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
        transpose_partial3((int16_t (*)[16])poly2, poly1                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
        transpose_partial3((int16_t (*)[16])poly2, poly1                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
        transpose_partial3((int16_t (*)[16])poly2, poly1                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                                   poly1 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("transpose partial 3: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        transpose_partial3_inv(poly2                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])poly1);
        transpose_partial3_inv(poly2                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])poly1);
        transpose_partial3_inv(poly2                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])poly1);
        transpose_partial3_inv(poly2                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])poly1);
        transpose_partial3_inv(poly2                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])poly1);
        transpose_partial3_inv(poly2                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])poly1);
        transpose_partial3_inv(poly2                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])poly1);
        transpose_partial3_inv(poly2                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                               poly2                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                               poly2 + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])poly1);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("transpose partial 3 inv: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);



}

