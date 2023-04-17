#ifndef __AVX2_BASEMUL_KARATSUBA_H
#define __AVX2_BASEMUL_KARATSUBA_H

#include <stdint.h>
#include <stddef.h>

#include "__avx2_basemul.h"

extern void __asm_karatsuba8(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);
extern void __asm_negacyclic_karatsuba8(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);

extern void __asm_weighted_karatsuba16(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff, int16_t *twiddle);
extern void __asm_weighted_double_karatsuba16(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff, int16_t *twiddle);

void cyclic_karatsuba2_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);
void negacyclic_karatsuba2_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);

void cyclic_karatsuba4_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);
void negacyclic_karatsuba4_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);

void cyclic_karatsuba8_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);
void negacyclic_karatsuba8_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2);

void karatsuba_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, size_t len);
void weighted_Karatsuba_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, size_t len, int16_t *twiddle);


static inline
void karatsuba2_Rmod_int16x16(int16x16_t des[3], int16x16_t src1[2], int16x16_t src2[2]){

    int16x16_t src1mid[1], src2mid[1];
    int16x16_t desmid[1];

    src1mid[0] = add_int16x16(src1[0], src1[1]);
    src2mid[0] = add_int16x16(src2[0], src2[1]);

    des[0] = montmulmod_int16x16(src1[0], src2[0]);
    des[2] = montmulmod_int16x16(src1[1], src2[1]);
    desmid[0] = montmulmod_int16x16(src1mid[0], src2mid[0]);

    desmid[0] = sub_int16x16(desmid[0], des[0]);
    desmid[0] = sub_int16x16(desmid[0], des[2]);

    des[1] = desmid[0];

}

static inline
void cyclic_karatsuba2_Rmod_int16x16(int16x16_t des[2], int16x16_t src1[2], int16x16_t src2[2]){

    int16x16_t src1mid[1], src2mid[1];
    int16x16_t desmid[1];
    int16x16_t deshi[1];

    src1mid[0] = add_int16x16(src1[0], src1[1]);
    src2mid[0] = add_int16x16(src2[0], src2[1]);

    desmid[0] = montmulmod_int16x16(src1mid[0], src2mid[0]);
    des[0] = montmulmod_int16x16(src1[0], src2[0]);
    deshi[0] = montmulmod_int16x16(src1[1], src2[1]);

    des[0] = add_int16x16(des[0], deshi[0]);
    des[1] = sub_int16x16(desmid[0], des[0]);

}

static inline
void negacyclic_karatsuba2_Rmod_int16x16(int16x16_t des[2], int16x16_t src1[2], int16x16_t src2[2]){

    int16x16_t src1mid[1], src2mid[1];
    int16x16_t desmid[1];
    int16x16_t deshi[1];

    src1mid[0] = add_int16x16(src1[0], src1[1]);
    src2mid[0] = add_int16x16(src2[0], src2[1]);

    desmid[0] = montmulmod_int16x16(src1mid[0], src2mid[0]);
    des[0] = montmulmod_int16x16(src1[0], src2[0]);
    deshi[0] = montmulmod_int16x16(src1[1], src2[1]);

    desmid[0] = sub_int16x16(desmid[0], des[0]);
    desmid[0] = sub_int16x16(desmid[0], deshi[0]);
    des[0] = sub_int16x16(des[0], deshi[0]);

    des[1] = desmid[0];

}

static inline
void karatsuba4_Rmod_int16x16(int16x16_t des[7], int16x16_t src1[4], int16x16_t src2[4]){

    int16x16_t src1mid[2], src2mid[2];
    int16x16_t desmid[3];

    src1mid[0] = add_int16x16(src1[0], src1[2]);
    src1mid[1] = add_int16x16(src1[1], src1[3]);

    src2mid[0] = add_int16x16(src2[0], src2[2]);
    src2mid[1] = add_int16x16(src2[1], src2[3]);

    schoolbook2_Rmod_int16x16(desmid, src1mid, src2mid);
    schoolbook2_Rmod_int16x16(des, src1, src2);
    schoolbook2_Rmod_int16x16(des + 4, src1 + 2, src2 + 2);

    desmid[0] = sub_int16x16(desmid[0], des[0]);
    desmid[0] = sub_int16x16(desmid[0], des[4]);
    desmid[1] = sub_int16x16(desmid[1], des[1]);
    desmid[1] = sub_int16x16(desmid[1], des[5]);
    desmid[2] = sub_int16x16(desmid[2], des[2]);
    desmid[2] = sub_int16x16(desmid[2], des[6]);

    des[2] = add_int16x16(des[2], desmid[0]);
    des[4] = add_int16x16(des[4], desmid[2]);
    des[3] = desmid[1];

}

static inline
void cyclic_karatsuba4_Rmod_int16x16(int16x16_t des[4], int16x16_t src1[4], int16x16_t src2[4]){

    int16x16_t src1mid[2], src2mid[2];
    int16x16_t desmid[3];
    int16x16_t deshi[3];

    src1mid[0] = add_int16x16(src1[0], src1[2]);
    src1mid[1] = add_int16x16(src1[1], src1[3]);

    src2mid[0] = add_int16x16(src2[0], src2[2]);
    src2mid[1] = add_int16x16(src2[1], src2[3]);

    schoolbook2_Rmod_int16x16(desmid, src1mid, src2mid);
    schoolbook2_Rmod_int16x16(des, src1, src2);
    schoolbook2_Rmod_int16x16(deshi, src1 + 2, src2 + 2);

    desmid[1] = sub_int16x16(desmid[1], des[1]);
    desmid[1] = sub_int16x16(desmid[1], deshi[1]);

    des[1] = add_int16x16(des[1], deshi[1]);

    des[3] = desmid[1];


    deshi[0] = add_int16x16(des[0], deshi[0]);
    deshi[2] = add_int16x16(des[2], deshi[2]);

    deshi[0] = barrett_int16x16(sub_int16x16(deshi[0], deshi[2]));

    des[0] = add_int16x16(desmid[2], deshi[0]);
    des[2] = sub_int16x16(desmid[0], deshi[0]);

}

static inline
void negacyclic_karatsuba4_Rmod_int16x16(int16x16_t des[4], int16x16_t src1[4], int16x16_t src2[4]){

    int16x16_t src1mid[2], src2mid[2];
    int16x16_t desmid[3];
    int16x16_t deshi[3];

    src1mid[0] = add_int16x16(src1[0], src1[2]);
    src1mid[1] = add_int16x16(src1[1], src1[3]);

    src2mid[0] = add_int16x16(src2[0], src2[2]);
    src2mid[1] = add_int16x16(src2[1], src2[3]);

    des[0] = montmulmod_int16x16(src1[0], src2[0]);
    des[2] = montmulmod_int16x16(src1[1], src2[1]);

    deshi[0] = montmulmod_int16x16(src1[2], src2[2]);
    deshi[2] = montmulmod_int16x16(src1[3], src2[3]);

    desmid[0] = montmulmod_int16x16(src1mid[0], src2mid[0]);
    desmid[2] = montmulmod_int16x16(src1mid[1], src2mid[1]);

    deshi[0] = sub_int16x16(des[2], deshi[0]);
    deshi[2] = add_int16x16(des[0], deshi[2]);


    desmid[1] = add_int16x16(montmulmod_int16x16(src1mid[0], src2mid[1]),
                             montmulmod_int16x16(src1mid[1], src2mid[0]));

    desmid[2] = sub_int16x16(deshi[0], desmid[2]);
    desmid[2] = add_int16x16(deshi[2], desmid[2]);


    des[1] = add_int16x16(montmulmod_int16x16(src1[0], src2[1]),
                          montmulmod_int16x16(src1[1], src2[0]));

    desmid[0] = add_int16x16(deshi[0], desmid[0]);
    desmid[0] = sub_int16x16(desmid[0], deshi[2]);

    deshi[1] = add_int16x16(montmulmod_int16x16(src1[2], src2[3]),
                            montmulmod_int16x16(src1[3], src2[2]));


    des[0] = barrett_int16x16(desmid[2]);
    des[2] = barrett_int16x16(desmid[0]);


    desmid[1] = sub_int16x16(desmid[1], des[1]);
    desmid[1] = sub_int16x16(desmid[1], deshi[1]);

    des[1] = sub_int16x16(des[1], deshi[1]);
    des[3] = desmid[1];

}

static inline
void negacyclic_karatsuba4_precompute_Rmod_int16x16(
    int16x16_t des[4], int16x16_t src1[4], int16x16_t src2[4], int16x16_t src2inv[4]){

    int16x16_t src1mid[2], src2mid[2];
    int16x16_t desmid[3];
    int16x16_t deshi[3];

    src1mid[0] = add_int16x16(src1[0], src1[2]);
    src1mid[1] = add_int16x16(src1[1], src1[3]);

    src2mid[0] = add_int16x16(src2[0], src2[2]);
    src2mid[1] = add_int16x16(src2[1], src2[3]);

    des[0] = montmulmod_precompute_int16x16(src1[0], src2[0], src2inv[0]);
    des[2] = montmulmod_precompute_int16x16(src1[1], src2[1], src2inv[1]);

    deshi[0] = montmulmod_precompute_int16x16(src1[2], src2[2], src2inv[2]);
    deshi[2] = montmulmod_precompute_int16x16(src1[3], src2[3], src2inv[3]);

    desmid[0] = montmulmod_int16x16(src1mid[0], src2mid[0]);
    desmid[2] = montmulmod_int16x16(src1mid[1], src2mid[1]);

    deshi[0] = sub_int16x16(des[2], deshi[0]);
    deshi[2] = add_int16x16(des[0], deshi[2]);


    desmid[1] = add_int16x16(montmulmod_int16x16(src1mid[0], src2mid[1]),
                             montmulmod_int16x16(src1mid[1], src2mid[0]));

    desmid[2] = sub_int16x16(deshi[0], desmid[2]);
    desmid[2] = add_int16x16(deshi[2], desmid[2]);


    des[1] = add_int16x16(montmulmod_precompute_int16x16(src1[0], src2[1], src2inv[1]),
                          montmulmod_precompute_int16x16(src1[1], src2[0], src2inv[0]));

    desmid[0] = add_int16x16(deshi[0], desmid[0]);
    desmid[0] = sub_int16x16(desmid[0], deshi[2]);

    deshi[1] = add_int16x16(montmulmod_precompute_int16x16(src1[2], src2[3], src2inv[3]),
                            montmulmod_precompute_int16x16(src1[3], src2[2], src2inv[2]));


    des[0] = barrett_int16x16(desmid[2]);
    des[2] = barrett_int16x16(desmid[0]);


    desmid[1] = sub_int16x16(desmid[1], des[1]);
    desmid[1] = sub_int16x16(desmid[1], deshi[1]);

    des[1] = sub_int16x16(des[1], deshi[1]);
    des[3] = desmid[1];

}

static inline
void karatsuba8_Rmod_int16x16(int16x16_t des[15], int16x16_t src1[8], int16x16_t src2[8]){

    int16x16_t src1mid[4], src2mid[4];
    int16x16_t desmid[7];

    for(size_t i = 0; i < 4; i++){
        src1mid[i] = add_int16x16(src1[i], src1[i + 4]);
    }
    for(size_t i = 0; i < 4; i++){
        src2mid[i] = add_int16x16(src2[i], src2[i + 4]);
    }

    schoolbook4_Rmod_int16x16(desmid, src1mid, src2mid);
    schoolbook4_Rmod_int16x16(des, src1, src2);
    for(size_t i = 0; i < 7; i++){
        desmid[i] = barrett_int16x16(sub_int16x16(desmid[i], des[i]));
    }
    schoolbook4_Rmod_int16x16(des + 8, src1 + 4, src2 + 4);
    for(size_t i = 0; i < 7; i++){
        desmid[i] = barrett_int16x16(sub_int16x16(desmid[i], des[i + 8]));
    }
    des[7] = desmid[3];
    for(size_t i = 0; i < 3; i++){
        des[4 + i] = barrett_int16x16(add_int16x16(des[4 + i], desmid[i]));
    }
    for(size_t i = 4; i < 7; i++){
        des[4 + i] = barrett_int16x16(add_int16x16(des[4 + i], desmid[i]));
    }

}

static inline
void karatsuba8_barrett_Rmod_int16x16(int16x16_t des[15], int16x16_t src1[8], int16x16_t src2[8]){

    int16x16_t src1mid[4], src2mid[4];
    int16x16_t desmid[7];


    for(size_t i = 0; i < 4; i++){
        src1mid[i] = barrett_int16x16(add_int16x16(src1[i], src1[i + 4]));
    }
    for(size_t i = 0; i < 4; i++){
        src2mid[i] = barrett_int16x16(add_int16x16(src2[i], src2[i + 4]));
    }

    schoolbook4_Rmod_int16x16(desmid, src1mid, src2mid);
    schoolbook4_Rmod_int16x16(des, src1, src2);
    for(size_t i = 0; i < 7; i++){
        desmid[i] = sub_int16x16(desmid[i], des[i]);
    }
    schoolbook4_Rmod_int16x16(des + 8, src1 + 4, src2 + 4);
    for(size_t i = 0; i < 7; i++){
        desmid[i] = sub_int16x16(desmid[i], des[i + 8]);
    }
    des[7] = barrett_int16x16(desmid[3]);
    for(size_t i = 0; i < 3; i++){
        des[4 + i] = barrett_int16x16(add_int16x16(des[4 + i], desmid[i]));
    }
    for(size_t i = 4; i < 7; i++){
        des[4 + i] = barrett_int16x16(add_int16x16(des[4 + i], desmid[i]));
    }

}

static inline
void cyclic_karatsuba8_Rmod_int16x16(int16x16_t des[8], int16x16_t src1[8], int16x16_t src2[8]){

    int16x16_t src1mid[4], src2mid[4];
    int16x16_t desmid[7];
    int16x16_t deshi[7];

    for(size_t i = 0; i < 4; i++){
        src1mid[i] = add_int16x16(src1[i], src1[i + 4]);
    }
    for(size_t i = 0; i < 4; i++){
        src2mid[i] = add_int16x16(src2[i], src2[i + 4]);
    }

    schoolbook4_Rmod_int16x16(desmid, src1mid, src2mid);
    schoolbook4_Rmod_int16x16(des, src1, src2);
    for(size_t i = 0; i < 7; i++){
        desmid[i] = sub_int16x16(desmid[i], des[i]);
    }
    schoolbook4_Rmod_int16x16(deshi, src1 + 4, src2 + 4);
    for(size_t i = 0; i < 7; i++){
        desmid[i] = sub_int16x16(desmid[i], deshi[i]);
    }

    des[7] = barrett_int16x16(desmid[3]);
    for(size_t i = 0; i < 3; i++){
        des[i] = barrett_int16x16(add_int16x16(add_int16x16(des[i], deshi[i]), desmid[i + 4]));
    }
    des[3] = barrett_int16x16(add_int16x16(des[3], deshi[3]));
    for(size_t i = 4; i < 7; i++){
        des[i] = barrett_int16x16(add_int16x16(add_int16x16(des[i], deshi[i]), desmid[i - 4]));
    }

}

static inline
void negacyclic_karatsuba8_Rmod_int16x16(int16x16_t des[8], int16x16_t src1[8], int16x16_t src2[8]){

    int16x16_t src1mid[4], src2mid[4];
    int16x16_t desmid[7];
    int16x16_t deshi[7];

    for(size_t i = 0; i < 4; i++){
        src1mid[i] = add_int16x16(src1[i], src1[i + 4]);
    }
    for(size_t i = 0; i < 4; i++){
        src2mid[i] = add_int16x16(src2[i], src2[i + 4]);
    }

    schoolbook4_Rmod_int16x16(desmid, src1mid, src2mid);
    schoolbook4_Rmod_int16x16(des, src1, src2);
    schoolbook4_Rmod_int16x16(deshi, src1 + 4, src2 + 4);

    des[7] = barrett_int16x16(sub_int16x16(sub_int16x16(desmid[3], des[3]), deshi[3]));
    des[3] = barrett_int16x16(sub_int16x16(des[3], deshi[3]));

    for(size_t i = 0; i < 3; i++){
        deshi[i] = sub_int16x16(des[i + 4], deshi[i]);
        deshi[i + 4] = add_int16x16(des[i], deshi[i + 4]);
        des[i]     = barrett_int16x16(
                            add_int16x16(sub_int16x16(deshi[i], desmid[i + 4]), deshi[i + 4]));
        des[i + 4] = barrett_int16x16(
                            sub_int16x16(add_int16x16(deshi[i], desmid[i]), deshi[i + 4]));
    }

}

static inline
void negacyclic_karatsuba8_precompute_Rmod_int16x16(
    int16x16_t des[8], int16x16_t src1[8], int16x16_t src2[8], int16x16_t src2inv[8]){

    int16x16_t src1mid[4], src2mid[4];
    int16x16_t desmid[7];
    int16x16_t deshi[7];

    for(size_t i = 0; i < 4; i++){
        src1mid[i] = add_int16x16(src1[i], src1[i + 4]);
    }
    for(size_t i = 0; i < 4; i++){
        src2mid[i] = add_int16x16(src2[i], src2[i + 4]);
    }

    schoolbook4_Rmod_int16x16(desmid, src1mid, src2mid);
    schoolbook4_precompute_Rmod_int16x16(des, src1, src2, src2inv);
    schoolbook4_precompute_Rmod_int16x16(deshi, src1 + 4, src2 + 4, src2inv + 4);

    des[7] = barrett_int16x16(sub_int16x16(sub_int16x16(desmid[3], des[3]), deshi[3]));
    des[3] = barrett_int16x16(sub_int16x16(des[3], deshi[3]));

    for(size_t i = 0; i < 3; i++){
        deshi[i] = sub_int16x16(des[i + 4], deshi[i]);
    }
    for(size_t i = 0; i < 3; i++){
        deshi[i + 4] = add_int16x16(des[i], deshi[i + 4]);
    }

    for(size_t i = 0; i < 3; i++){
        des[i]     = barrett_int16x16(
                            add_int16x16(sub_int16x16(deshi[i], desmid[i + 4]), deshi[i + 4]));
        des[i + 4] = barrett_int16x16(
                            sub_int16x16(add_int16x16(deshi[i], desmid[i]), deshi[i + 4]));
    }

}

static inline
void karatsuba16_Rmod_int16x16(int16x16_t des[31], int16x16_t src1[16], int16x16_t src2[16]){

    int16x16_t src1mid[8], src2mid[8];
    int16x16_t desmid[15];

    for(size_t i = 0; i < 8; i++){
        src1mid[i] = add_int16x16(src1[i], src1[i + 8]);
    }

    for(size_t i = 0; i < 8; i++){
        src2mid[i] = add_int16x16(src2[i], src2[i + 8]);
    }

    karatsuba8_barrett_Rmod_int16x16(desmid, src1mid, src2mid);
    karatsuba8_Rmod_int16x16(des, src1, src2);
    for(size_t i = 0; i < 15; i++){
        desmid[i] = sub_int16x16(desmid[i], des[i]);
    }
    karatsuba8_Rmod_int16x16(des + 16, src1 + 8, src2 + 8);
    for(size_t i = 0; i < 15; i++){
        desmid[i] = sub_int16x16(desmid[i], des[i + 16]);
    }
    des[15] = barrett_int16x16(desmid[7]);
    for(size_t i = 0; i < 7; i++){
        des[8 + i] = barrett_int16x16(add_int16x16(des[8 + i], desmid[i]));
    }
    for(size_t i = 8; i < 15; i++){
        des[8 + i] = barrett_int16x16(add_int16x16(des[8 + i], desmid[i]));
    }

}



#endif

