#ifndef __AVX2_BASEMUL_FFT_H
#define __AVX2_BASEMUL_FFT_H

#include <stdint.h>

#include "__avx2_wrap.h"
#include "__avx2_basemul.h"
#include "__avx2_basemul_Karatsuba.h"

extern void __asm_cyclic_FFT2(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);
extern void __asm_negacyclic_FFT4(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);

extern void __asm_cyclic_FFT4(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);
extern void __asm_cyclic_FFT8(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);

extern void __asm_cyclic_FFT16_precompute(int16_t *des, int16_t *src1, int16_t *src2, int16_t *src2inv, int16_t *_const_buff);

void cyclic_FFT2_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);
void cyclic_FFT4_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);
void cyclic_FFT8_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);

void cyclic_FFT16_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);
void cyclic_FFT16_precompute_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, int16_t *src2inv);

void weighted16_Rmod_FFT2_int16x16(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_CT_table, int16_t *_GS_table);

static inline
void cyclic_FFT2_Rmod_int16x16(int16x16_t des[2], int16x16_t src1[2], int16x16_t src2[2]){

    int16x16_t t1[2], t2[2];

    t1[0] = add_int16x16(src1[0], src1[1]);
    t1[1] = sub_int16x16(src1[0], src1[1]);

    t2[0] = add_int16x16(src2[0], src2[1]);
    t2[1] = sub_int16x16(src2[0], src2[1]);

    t1[0] = montmulmod_int16x16(t1[0], t2[0]);
    t1[1] = montmulmod_int16x16(t1[1], t2[1]);

    des[0] = add_int16x16(t1[0], t1[1]);
    des[1] = sub_int16x16(t1[0], t1[1]);

}

static inline
void cyclic_FFT4_Rmod_int16x16(int16x16_t des[4], int16x16_t src1[4], int16x16_t src2[4]){

    int16x16_t t1[4], t2[4];
    int16x16_t buff[4];

    t1[0] = add_int16x16(src1[0], src1[2]);
    t1[2] = sub_int16x16(src1[0], src1[2]);
    t1[1] = add_int16x16(src1[1], src1[3]);
    t1[3] = sub_int16x16(src1[1], src1[3]);


    t2[0] = add_int16x16(src2[0], src2[2]);
    t2[2] = sub_int16x16(src2[0], src2[2]);
    t2[1] = add_int16x16(src2[1], src2[3]);
    t2[3] = sub_int16x16(src2[1], src2[3]);

    buff[0] = add_int16x16(t1[0], t1[1]);
    buff[1] = sub_int16x16(t1[0], t1[1]);
    buff[2] = add_int16x16(t2[0], t2[1]);
    buff[3] = sub_int16x16(t2[0], t2[1]);

    buff[0] = montmulmod_int16x16(buff[0], buff[2]);
    buff[1] = montmulmod_int16x16(buff[1], buff[3]);

    t1[0] = add_int16x16(buff[0], buff[1]);
    t1[1] = sub_int16x16(buff[0], buff[1]);

    negacyclic_schoolbook2_Rmod_int16x16(buff + 2, t1 + 2, t2 + 2);

    t1[2] = add_int16x16(buff[2], buff[2]);
    t1[3] = add_int16x16(buff[3], buff[3]);

    des[0] = barrett_int16x16(add_int16x16(t1[0], t1[2]));
    des[2] = barrett_int16x16(sub_int16x16(t1[0], t1[2]));
    des[1] = barrett_int16x16(add_int16x16(t1[1], t1[3]));
    des[3] = barrett_int16x16(sub_int16x16(t1[1], t1[3]));

}

static inline
void cyclic_FFT8_Rmod_int16x16(int16x16_t des[8], int16x16_t src1[8], int16x16_t src2[8]){

    int16x16_t t1[8], t2[8];
    int16x16_t buff[8];

    t1[0] = add_int16x16(src1[0], src1[4]);
    t1[4] = sub_int16x16(src1[0], src1[4]);
    t1[1] = add_int16x16(src1[1], src1[5]);
    t1[5] = sub_int16x16(src1[1], src1[5]);
    t1[2] = add_int16x16(src1[2], src1[6]);
    t1[6] = sub_int16x16(src1[2], src1[6]);
    t1[3] = add_int16x16(src1[3], src1[7]);
    t1[7] = sub_int16x16(src1[3], src1[7]);

    src1[0] = add_int16x16(t1[0], t1[2]);
    src1[2] = sub_int16x16(t1[0], t1[2]);

    src1[1] = add_int16x16(t1[1], t1[3]);
    src1[3] = sub_int16x16(t1[1], t1[3]);
    src1[1] = barrett_int16x16(src1[1]);

    t1[0] = add_int16x16(src1[0], src1[1]);
    t1[1] = sub_int16x16(src1[0], src1[1]);

    t2[0] = add_int16x16(src2[0], src2[4]);
    t2[4] = sub_int16x16(src2[0], src2[4]);
    t2[1] = add_int16x16(src2[1], src2[5]);
    t2[5] = sub_int16x16(src2[1], src2[5]);
    t2[2] = add_int16x16(src2[2], src2[6]);
    t2[6] = sub_int16x16(src2[2], src2[6]);
    t2[3] = add_int16x16(src2[3], src2[7]);
    t2[7] = sub_int16x16(src2[3], src2[7]);

    src2[0] = add_int16x16(t2[0], t2[2]);
    src2[2] = sub_int16x16(t2[0], t2[2]);

    src2[1] = add_int16x16(t2[1], t2[3]);
    src2[3] = sub_int16x16(t2[1], t2[3]);
    src2[1] = barrett_int16x16(src2[1]);

    t2[0] = add_int16x16(src2[0], src2[1]);
    t2[1] = sub_int16x16(src2[0], src2[1]);

    buff[0] = montmulmod_int16x16(t1[0], t2[0]);
    buff[1] = montmulmod_int16x16(t1[1], t2[1]);
    negacyclic_schoolbook2_Rmod_int16x16(buff + 2, src1 + 2, src2 + 2);
    buff[2] = add_int16x16(buff[2], buff[2]);
    buff[3] = add_int16x16(buff[3], buff[3]);

    t1[0] = add_int16x16(buff[0], buff[1]);
    t1[1] = sub_int16x16(buff[0], buff[1]);

    buff[0] = add_int16x16(t1[0], buff[2]);
    buff[2] = sub_int16x16(t1[0], buff[2]);
    buff[1] = add_int16x16(t1[1], buff[3]);
    buff[3] = sub_int16x16(t1[1], buff[3]);

    negacyclic_karatsuba4_Rmod_int16x16(buff + 4, t1 + 4, t2 + 4);
    buff[4] = add_int16x16(buff[4], buff[4]);
    buff[5] = add_int16x16(buff[5], buff[5]);
    buff[6] = add_int16x16(buff[6], buff[6]);
    buff[7] = add_int16x16(buff[7], buff[7]);

    buff[4] = add_int16x16(buff[4], buff[4]);
    buff[5] = add_int16x16(buff[5], buff[5]);
    buff[6] = add_int16x16(buff[6], buff[6]);
    buff[7] = add_int16x16(buff[7], buff[7]);

    buff[4] = barrett_int16x16(buff[4]);
    buff[5] = barrett_int16x16(buff[5]);
    buff[6] = barrett_int16x16(buff[6]);
    buff[7] = barrett_int16x16(buff[7]);

    des[0] = add_int16x16(buff[0], buff[4]);
    des[4] = sub_int16x16(buff[0], buff[4]);
    des[1] = add_int16x16(buff[1], buff[5]);
    des[5] = sub_int16x16(buff[1], buff[5]);
    des[2] = add_int16x16(buff[2], buff[6]);
    des[6] = sub_int16x16(buff[2], buff[6]);
    des[3] = add_int16x16(buff[3], buff[7]);
    des[7] = sub_int16x16(buff[3], buff[7]);


}

// TODO: don't modify inputs
static inline
void cyclic_FFT16_Rmod_int16x16(int16x16_t des[16], int16x16_t src1[16], int16x16_t src2[16]){

    int16x16_t t1[16], t2[16];
    int16x16_t buff[16];

    for(size_t i = 0; i < 2; i++){

        t1[i +  0] = add_int16x16(src1[i +  0], src1[i +  8]);
        t1[i +  8] = sub_int16x16(src1[i +  0], src1[i +  8]);
        t1[i +  2] = add_int16x16(src1[i +  2], src1[i + 10]);
        t1[i + 10] = sub_int16x16(src1[i +  2], src1[i + 10]);
        t1[i +  4] = add_int16x16(src1[i +  4], src1[i + 12]);
        t1[i + 12] = sub_int16x16(src1[i +  4], src1[i + 12]);
        t1[i +  6] = add_int16x16(src1[i +  6], src1[i + 14]);
        t1[i + 14] = sub_int16x16(src1[i +  6], src1[i + 14]);

        src1[i +  0] = add_int16x16(t1[i +  0], t1[i +  4]);
        src1[i +  4] = sub_int16x16(t1[i +  0], t1[i +  4]);
        src1[i +  2] = add_int16x16(t1[i +  2], t1[i +  6]);
        src1[i +  6] = sub_int16x16(t1[i +  2], t1[i +  6]);

        src1[i +  0] = barrett_int16x16(src1[i +  0]);
        src1[i +  2] = barrett_int16x16(src1[i +  2]);
        src1[i +  4] = barrett_int16x16(src1[i +  4]);
        src1[i +  6] = barrett_int16x16(src1[i +  6]);

        t1[i +  0] = add_int16x16(src1[i +  0], src1[i +  2]);
        t1[i +  2] = sub_int16x16(src1[i +  0], src1[i +  2]);
    }

    for(size_t i = 0; i < 2; i++){

        t2[i +  0] = add_int16x16(src2[i +  0], src2[i +  8]);
        t2[i +  8] = sub_int16x16(src2[i +  0], src2[i +  8]);
        t2[i +  2] = add_int16x16(src2[i +  2], src2[i + 10]);
        t2[i + 10] = sub_int16x16(src2[i +  2], src2[i + 10]);
        t2[i +  4] = add_int16x16(src2[i +  4], src2[i + 12]);
        t2[i + 12] = sub_int16x16(src2[i +  4], src2[i + 12]);
        t2[i +  6] = add_int16x16(src2[i +  6], src2[i + 14]);
        t2[i + 14] = sub_int16x16(src2[i +  6], src2[i + 14]);

        src2[i +  0] = add_int16x16(t2[i +  0], t2[i +  4]);
        src2[i +  4] = sub_int16x16(t2[i +  0], t2[i +  4]);
        src2[i +  2] = add_int16x16(t2[i +  2], t2[i +  6]);
        src2[i +  6] = sub_int16x16(t2[i +  2], t2[i +  6]);

        src2[i +  0] = barrett_int16x16(src2[i +  0]);
        src2[i +  2] = barrett_int16x16(src2[i +  2]);
        src2[i +  4] = barrett_int16x16(src2[i +  4]);
        src2[i +  6] = barrett_int16x16(src2[i +  6]);

        t2[i +  0] = add_int16x16(src2[i +  0], src2[i +  2]);
        t2[i +  2] = sub_int16x16(src2[i +  0], src2[i +  2]);
    }

    src1[0] = add_int16x16(t1[0], t1[1]);
    src1[1] = sub_int16x16(t1[0], t1[1]);
    src2[0] = add_int16x16(t2[0], t2[1]);
    src2[1] = sub_int16x16(t2[0], t2[1]);

    t1[0] = montmulmod_int16x16(src1[0], src2[0]);
    t1[1] = montmulmod_int16x16(src1[1], src2[1]);

    buff[0] = add_int16x16(t1[0], t1[1]);
    buff[1] = sub_int16x16(t1[0], t1[1]);

    negacyclic_schoolbook2_Rmod_int16x16(buff + 2, t1 + 2, t2 + 2);

    buff[2] = add_int16x16(buff[2], buff[2]);
    buff[3] = add_int16x16(buff[3], buff[3]);

    t1[0] = add_int16x16(buff[0], buff[2]);
    t1[2] = sub_int16x16(buff[0], buff[2]);
    t1[1] = add_int16x16(buff[1], buff[3]);
    t1[3] = sub_int16x16(buff[1], buff[3]);

    negacyclic_karatsuba4_Rmod_int16x16(t1 + 4, src1 + 4, src2 + 4);

    for(size_t i = 4; i < 8; i++){
        t1[i] = barrett_int16x16(add_int16x16(t1[i], t1[i]));
        t1[i] = barrett_int16x16(add_int16x16(t1[i], t1[i]));
    }

    for(size_t i = 0; i < 4; i++){
        buff[i + 4] = sub_int16x16(t1[i], t1[i + 4]);
        buff[i    ] = add_int16x16(t1[i], t1[i + 4]);
    }

    negacyclic_karatsuba8_Rmod_int16x16(buff + 8, t1 + 8, t2 + 8);

    for(size_t i = 8; i < 16; i++){
        buff[i] = barrett_int16x16(add_int16x16(buff[i], buff[i]));
        buff[i] = add_int16x16(buff[i], buff[i]);
        buff[i] = barrett_int16x16(add_int16x16(buff[i], buff[i]));
    }

    for(size_t i = 0; i < 8; i++){
        des[i    ] = barrett_int16x16(add_int16x16(buff[i], buff[i + 8]));
        des[i + 8] = barrett_int16x16(sub_int16x16(buff[i], buff[i + 8]));
    }


}

// TODO: don't modify inputs
static inline
void cyclic_FFT16_precompute_int16x16(int16x16_t des[16], int16x16_t src1[16], int16x16_t src2[16], int16x16_t src2inv[16]){

    int16x16_t t1[16];
    int16x16_t buff[16];

    for(size_t i = 0; i < 2; i++){

        t1[i +  0] = add_int16x16(src1[i +  0], src1[i +  8]);
        t1[i +  8] = sub_int16x16(src1[i +  0], src1[i +  8]);
        t1[i +  2] = add_int16x16(src1[i +  2], src1[i + 10]);
        t1[i + 10] = sub_int16x16(src1[i +  2], src1[i + 10]);
        t1[i +  4] = add_int16x16(src1[i +  4], src1[i + 12]);
        t1[i + 12] = sub_int16x16(src1[i +  4], src1[i + 12]);
        t1[i +  6] = add_int16x16(src1[i +  6], src1[i + 14]);
        t1[i + 14] = sub_int16x16(src1[i +  6], src1[i + 14]);

        src1[i +  0] = add_int16x16(t1[i +  0], t1[i +  4]);
        src1[i +  4] = sub_int16x16(t1[i +  0], t1[i +  4]);
        src1[i +  2] = add_int16x16(t1[i +  2], t1[i +  6]);
        src1[i +  6] = sub_int16x16(t1[i +  2], t1[i +  6]);

        src1[i +  0] = barrett_int16x16(src1[i +  0]);
        src1[i +  2] = barrett_int16x16(src1[i +  2]);
        src1[i +  4] = barrett_int16x16(src1[i +  4]);
        src1[i +  6] = barrett_int16x16(src1[i +  6]);

        t1[i +  0] = add_int16x16(src1[i +  0], src1[i +  2]);
        t1[i +  2] = sub_int16x16(src1[i +  0], src1[i +  2]);

    }

    src1[0] = add_int16x16(t1[0], t1[1]);
    src1[1] = sub_int16x16(t1[0], t1[1]);

    t1[0] = montmulmod_precompute_int16x16(src1[0], src2[0], src2inv[0]);
    t1[1] = montmulmod_precompute_int16x16(src1[1], src2[1], src2inv[1]);

    buff[0] = add_int16x16(t1[0], t1[1]);
    buff[1] = sub_int16x16(t1[0], t1[1]);

    negacyclic_schoolbook2_precompute_Rmod_int16x16(buff + 2, t1 + 2, src2 + 2, src2inv + 2);

    t1[0] = add_int16x16(buff[0], buff[2]);
    t1[2] = sub_int16x16(buff[0], buff[2]);
    t1[1] = add_int16x16(buff[1], buff[3]);
    t1[3] = sub_int16x16(buff[1], buff[3]);

    negacyclic_karatsuba4_precompute_Rmod_int16x16(t1 + 4, src1 + 4, src2 + 4, src2inv + 4);

    for(size_t i = 0; i < 4; i++){
        buff[i + 4] = sub_int16x16(t1[i], t1[i + 4]);
        buff[i    ] = add_int16x16(t1[i], t1[i + 4]);
    }

    negacyclic_karatsuba8_precompute_Rmod_int16x16(buff + 8, t1 + 8, src2 + 8, src2inv + 8);

    for(size_t i = 0; i < 8; i++){
        des[i    ] = barrett_int16x16(add_int16x16(buff[i], buff[i + 8]));
        des[i + 8] = barrett_int16x16(sub_int16x16(buff[i], buff[i + 8]));
    }


}

#endif

