
#include "vec.h"
#include "ring.h"

void add_int16x16_(int16_t des[16], int16_t src1[16], int16_t src2[16]){

    for(size_t i = 0; i < 16; i++){
        coeff_ring.addZ(des + i, src1 + i, src2 + i);
    }

}

void sub_int16x16_(int16_t des[16], int16_t src1[16], int16_t src2[16]){

    for(size_t i = 0; i < 16; i++){
        coeff_ring.subZ(des + i, src1 + i, src2 + i);
    }

}

void mul_int16x16_(int16_t des[16], int16_t src1[16], int16_t src2[16]){

    for(size_t i = 0; i < 16; i++){
        coeff_ring.mulZ(des + i, src1 + i, src2 + i);
    }

}

void karatsuba_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, size_t len){

    if(len == 1){
        mul_int16x16_(des, src1, src2);
        return;
    }

    int16_t src1mid[(len / 2 ) * 16], src2mid[(len / 2 ) * 16];
    int16_t desmid[(len - 1 ) * 16];

    for(size_t i = 0; i < (len / 2); i++){
        add_int16x16_(src1mid + i * 16, src1 + (0 + i) * 16, src1 + ((len / 2) + i) * 16);
        add_int16x16_(src2mid + i * 16, src2 + (0 + i) * 16, src2 + ((len / 2) + i) * 16);
    }

    for(size_t i = 0; i < 2 * len - 1; i++){
        memset(des + i * 16, 0, 16 * sizeof(int16_t));
    }

    karatsuba_int16x16_(des, src1, src2, len / 2);
    karatsuba_int16x16_(des + len * 16, src1 + (len / 2) * 16, src2 + (len / 2) * 16, len / 2);
    karatsuba_int16x16_(desmid, src1mid, src2mid, len / 2);

    for(size_t i = 0; i < len - 1; i++){
        sub_int16x16_(desmid + i * 16, desmid + i * 16, des + i * 16);
        sub_int16x16_(desmid + i * 16, desmid + i * 16, des + (i + len) * 16);
    }

    for(size_t i = 0; i < len - 1; i++){
        add_int16x16_(des + (i + (len / 2)) * 16, des + (i + (len / 2)) * 16, desmid + i * 16);
    }

}



void weighted_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, size_t len, int16_t *twiddle){

    int16_t buff[(2 * len - 1) * 16 * sizeof (int16_t)];

    karatsuba_int16x16_(buff, src1, src2, len);

    for(size_t i = len; i < 2 * len - 1; i++){
        mul_int16x16_(buff + i * 16, buff + i * 16, twiddle);
        add_int16x16_(buff + (i - len) * 16, buff + (i - len) * 16, buff + i * 16);
    }

    memmove(des, buff, len * 16 * sizeof(int16_t));

}

