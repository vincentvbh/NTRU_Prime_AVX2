
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "params.h"
#include "NTT_params.h"

#include "tools.h"
#include "naive_mult.h"
#include "gen_table.h"
#include "ntt_c.h"

#include "ring.h"
#include "vec.h"
#include "__avx2.h"
#include "__avx2_wrap.h"
#include "__avx2_permute.h"
#include "rader.h"

int main(void){

    int16_t poly1[ARRAY_N], poly2[ARRAY_N];
    int16_t ref[ARRAY_N], res[ARRAY_N];
    int16_t poly1_NTT[7 * 256];
    int16_t poly2_NTT[7 * 256];
    int16_t res_NTT[7 * 256];
    int16_t poly2inv[ARRAY_N];

    int16_t buff1[256], buff2[256], buff3[256];

    int16_t scale, omega, twiddle, t, v;

    (void)omega;
    (void)res_NTT;
    (void)buff1;

    for(size_t i = 0; i < p; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

    for(size_t i = p; i < ARRAY_N; i++){
        poly1[i] = poly2[i] = 0;
    }

// ================================
// Barrett

    for(t = -32768; t < 0; t++){
        buff1[0] = t;
        __barrett_int16x16(buff1);
        v = buff1[0];
        if( (v < -BARRETT_BOUND) || (v > BARRETT_BOUND) ){
            printf("%8d: %8d\n", t, v);
        }
    }
    for(t = 32767; t > 0; t--){
        buff1[0] = t;
        __barrett_int16x16(buff1);
        v = buff1[0];
        if( (v < -BARRETT_BOUND) || (v > BARRETT_BOUND) ){
            printf("%8d: %8d\n", t, v);
        }
    }


// ================================
// 2 x 2

    printf("================================\n");

    for(size_t i = 2; i < ARRAY_N; i++){
        poly1_NTT[i] = poly1[i];
        poly2_NTT[i] = poly2[i];
    }
    for(size_t i = 2; i < ARRAY_N; i++){
        poly1[i] = 0;
        poly2[i] = 0;
    }

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 4, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_schoolbook2(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    for(size_t i = 2; i < ARRAY_N; i++){
        poly1[i] = poly1_NTT[i];
        poly2[i] = poly2_NTT[i];
    }

    scale = RmodQ;
    for(size_t i = 0; i < 3; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 3; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm schoolbook 2x2 passed!\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 2, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_cyclic_schoolbook2(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 2; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 2; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm schoolbook cyclic 2x2 passed!\n");

    twiddle = -1;
    naive_mulR(ref, poly1, poly2, 2, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_negacyclic_schoolbook2(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 2; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 2; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm schoolbook negacyclic 2x2 passed!\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 2, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_cyclic_FFT2(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 2; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    scale = 2;
    for(size_t i = 0; i < 2; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < 2; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm FFT cyclic 2x2 passed!\n");

#if 0

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 2, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    cyclic_FFT2_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 2; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    scale = 2;
    for(size_t i = 0; i < 2; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < 2; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("FFT cyclic 2x2 passed!\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 2, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    cyclic_karatsuba2_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 2; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 2; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("Karatsuba cyclic 2x2 passed!\n");

    twiddle = -1;
    naive_mulR(ref, poly1, poly2, 2, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    negacyclic_karatsuba2_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 2; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 2; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("Karatsuba negacyclic 2x2 passed!\n");

#endif


// ================================
// 4 x 4

    printf("================================\n");

    for(size_t i = 4; i < ARRAY_N; i++){
        poly1_NTT[i] = poly1[i];
        poly2_NTT[i] = poly2[i];
    }
    for(size_t i = 4; i < ARRAY_N; i++){
        poly1[i] = 0;
        poly2[i] = 0;
    }

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 8, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_schoolbook4(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    for(size_t i = 4; i < ARRAY_N; i++){
        poly1[i] = poly1_NTT[i];
        poly2[i] = poly2_NTT[i];
    }

    scale = RmodQ;
    for(size_t i = 0; i < 7; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 7; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm schoolbook 4x4 passed!\n");

    for(size_t i = 0; i < 256; i++){
        res[i] = 0;
    }

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 4, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_cyclic_schoolbook4(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 4; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm schoolbook cyclic 4x4 passed!\n");

    twiddle = -1;
    naive_mulR(ref, poly1, poly2, 4, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_negacyclic_schoolbook4(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 4; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm schoolbook negacyclic 4x4 passed!\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 4, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_cyclic_FFT4(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    scale = 2;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < 4; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm FFT cyclic 4x4 passed!\n");

    twiddle = -1;
    naive_mulR(ref, poly1, poly2, 4, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_negacyclic_FFT4(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    scale = 2;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < 4; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm FFT negacyclic 4x4 passed!\n");

#if 0

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 4, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    cyclic_karatsuba4_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 4; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("Karatsuba cyclic 4x4 passed!\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 4, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    cyclic_FFT4_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    scale = 4;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < 4; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("FFT cyclic 4x4 passed!\n");

    twiddle = -1;
    naive_mulR(ref, poly1, poly2, 4, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    negacyclic_karatsuba4_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 4; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 4; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("Karatsuba negacyclic 4x4 passed!\n");

#endif

// ================================
// 8 x 8

    printf("================================\n");

    for(size_t i = 8; i < ARRAY_N; i++){
        poly1_NTT[i] = poly1[i];
        poly2_NTT[i] = poly2[i];
    }
    for(size_t i = 8; i < ARRAY_N; i++){
        poly1[i] = 0;
        poly2[i] = 0;
    }

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 16, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_karatsuba8(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    for(size_t i = 8; i < ARRAY_N; i++){
        poly1[i] = poly1_NTT[i];
        poly2[i] = poly2_NTT[i];
    }

    scale = RmodQ;
    for(size_t i = 0; i < 15; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 15; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm Karatsuba 8x8 passed!\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 8, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_cyclic_FFT8(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 8; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    scale = 4;
    for(size_t i = 0; i < 8; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < 8; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm FFT cyclic 8x8 passed!\n");

    twiddle = -1;
    naive_mulR(ref, poly1, poly2, 8, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    __asm_negacyclic_karatsuba8(res, poly1, poly2, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 8; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 8; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm Karatsuba negacyclic 8x8 passed!\n");

#if 0

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 8, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    cyclic_karatsuba8_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 8; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 8; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("Karatsuba cyclic 8x8 passed!\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 8, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    cyclic_FFT8_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 8; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    scale = 8;
    for(size_t i = 0; i < 8; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < 8; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("FFT cyclic 8x8 passed!\n");

    twiddle = -1;
    naive_mulR(ref, poly1, poly2, 8, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    negacyclic_karatsuba8_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 8; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 8; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("Karatsuba negacyclic 8x8 passed!\n");

#endif

// ================================
// 16 x 16

    printf("================================\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 16, &twiddle, coeff_ring);

    memmove(buff2, poly2, 16 * sizeof(int16_t));
    memmove(buff3, poly2, 16 * sizeof(int16_t));

    for(size_t i = 0; i < 8; i++){
        coeff_ring.addZ(poly2 + i    , buff3 + i, buff3 + i + 8);
        coeff_ring.subZ(poly2 + i + 8, buff3 + i, buff3 + i + 8);
    }
    for(size_t i = 0; i < 4; i++){
        coeff_ring.addZ(buff3 + i    , poly2 + i, poly2 + i + 4);
        coeff_ring.subZ(buff3 + i + 4, poly2 + i, poly2 + i + 4);
    }
    for(size_t i = 0; i < 2; i++){
        coeff_ring.addZ(poly2 + i    , buff3 + i, buff3 + i + 2);
        coeff_ring.subZ(poly2 + i + 2, buff3 + i, buff3 + i + 2);
    }
    coeff_ring.addZ(buff3 + 0, poly2 + 0, poly2 + 1);
    coeff_ring.subZ(buff3 + 1, poly2 + 0, poly2 + 1);

    poly2[0] = buff3[0];
    poly2[1] = buff3[1];
    poly2[4] = buff3[4];
    poly2[5] = buff3[5];
    poly2[6] = buff3[6];
    poly2[7] = buff3[7];

    scale = 2;
    for(size_t i = 2; i < 4; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }
    scale = 4;
    for(size_t i = 4; i < 8; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }
    scale = 8;
    for(size_t i = 8; i < 16; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }

    scale = RmodQ;
    for(size_t i = 0; i < 16; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }

    scale = INV16;
    for(size_t i = 0; i < 16; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }

    scale = Qinv;
    for(size_t i = 0; i < 16; i++){
        poly2inv[i] = poly2[i] * scale;
    }

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    transpose((int16_t (*)[16])(poly2inv), (int16_t (*)[16])(poly2inv));
    __asm_cyclic_FFT16_precompute(res, poly1, poly2, poly2inv, const_buff);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    memmove(poly2, buff2, 16 * sizeof(int16_t));

    for(size_t i = 0; i < 16; i++){
        coeff_ring.memberZ(res + i, res + i);
    }

    for(size_t i = 0; i < 16; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("asm FFT cyclic precompute 16x16 passed!\n");

    transpose((int16_t (*)[16])(poly1),
              (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2),
              (int16_t (*)[16])(poly2));

    __asm_weighted_karatsuba16(ref, poly1, poly2, const_buff, const_twiddle51p);
    weighted16_Rmod_FFT2_int16x16(res, poly1, poly2, twiddle_CT, twiddle_GS);

    transpose((int16_t (*)[16])(ref),
              (int16_t (*)[16])(ref));

    transpose((int16_t (*)[16])(res),
              (int16_t (*)[16])(res));

    for(size_t i = 0; i < 16 * 16; i++){
        coeff_ring.memberZ(ref + i, ref + i);
    }

    scale = INV2;
    for(size_t i = 0; i < 16 * 16; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < 16; i++){
        for(size_t j = 0; j < 16; j++){
            if(ref[i * 16 + j] != res[i * 16 + j]){
                printf("%4zu, %4zu: %12d, %12d\n", i, j, ref[i * 16 + j], res[i * 16 + j]);
            }
        }
    }

    printf("asm weighted Karatsuba and FFT2 are compatible (16x16)!\n");

#if 0

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 16, &twiddle, coeff_ring);

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    cyclic_FFT16_Rmod_int16x16_(res, poly1, poly2);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    scale = RmodQ;
    for(size_t i = 0; i < 16; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    scale = 16;
    for(size_t i = 0; i < 16; i++){
        coeff_ring.mulZ(ref + i, ref + i, &scale);
    }

    for(size_t i = 0; i < 16; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("FFT cyclic 16x16 passed!\n");

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, 16, &twiddle, coeff_ring);

    memmove(buff2, poly2, 16 * sizeof(int16_t));

    memmove(buff3, poly2, 16 * sizeof(int16_t));

    for(size_t i = 0; i < 8; i++){
        coeff_ring.addZ(poly2 + i    , buff3 + i, buff3 + i + 8);
        coeff_ring.subZ(poly2 + i + 8, buff3 + i, buff3 + i + 8);
    }
    for(size_t i = 0; i < 4; i++){
        coeff_ring.addZ(buff3 + i    , poly2 + i, poly2 + i + 4);
        coeff_ring.subZ(buff3 + i + 4, poly2 + i, poly2 + i + 4);
    }
    for(size_t i = 0; i < 2; i++){
        coeff_ring.addZ(poly2 + i    , buff3 + i, buff3 + i + 2);
        coeff_ring.subZ(poly2 + i + 2, buff3 + i, buff3 + i + 2);
    }
    coeff_ring.addZ(buff3 + 0, poly2 + 0, poly2 + 1);
    coeff_ring.subZ(buff3 + 1, poly2 + 0, poly2 + 1);

    poly2[0] = buff3[0];
    poly2[1] = buff3[1];
    poly2[4] = buff3[4];
    poly2[5] = buff3[5];
    poly2[6] = buff3[6];
    poly2[7] = buff3[7];

    scale = 2;
    for(size_t i = 2; i < 4; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }
    scale = 4;
    for(size_t i = 4; i < 8; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }
    scale = 8;
    for(size_t i = 8; i < 16; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }

    scale = RmodQ;
    for(size_t i = 0; i < 16; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }

    scale = INV16;
    for(size_t i = 0; i < 16; i++){
        coeff_ring.mulZ(poly2 + i, poly2 + i, &scale);
    }

    scale = Qinv;
    for(size_t i = 0; i < 16; i++){
        poly2inv[i] = poly2[i] * scale;
    }

    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));
    transpose((int16_t (*)[16])(poly2inv), (int16_t (*)[16])(poly2inv));
    cyclic_FFT16_precompute_int16x16_(res, poly1, poly2, poly2inv);
    transpose((int16_t (*)[16])(res), (int16_t (*)[16])(res));
    transpose((int16_t (*)[16])(poly1), (int16_t (*)[16])(poly1));
    transpose((int16_t (*)[16])(poly2), (int16_t (*)[16])(poly2));

    memmove(poly2, buff2, 16 * sizeof(int16_t));

    for(size_t i = 0; i < 16; i++){
        coeff_ring.memberZ(res + i, res + i);
    }

    for(size_t i = 0; i < 16; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("FFT cyclic precompute 16x16 passed!\n");

#endif

// ================================

    printf("================================\n");


}

