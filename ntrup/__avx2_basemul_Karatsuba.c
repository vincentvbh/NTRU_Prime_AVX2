
#include <memory.h>

#include "__avx2_basemul.h"
#include "__avx2_basemul_Karatsuba.h"

#include "vec.h"

void cyclic_karatsuba2_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[2], src2x16[2];
    int16x16_t desx16[2];

    for(size_t i = 0; i < 2; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    cyclic_karatsuba2_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 2; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}

void negacyclic_karatsuba2_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[2], src2x16[2];
    int16x16_t desx16[2];

    for(size_t i = 0; i < 2; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    negacyclic_karatsuba2_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 2; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}

void cyclic_karatsuba4_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[4], src2x16[4];
    int16x16_t desx16[4];

    for(size_t i = 0; i < 4; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    cyclic_karatsuba4_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 4; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}

void negacyclic_karatsuba4_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[4], src2x16[4];
    int16x16_t desx16[4];

    for(size_t i = 0; i < 4; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    negacyclic_karatsuba4_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 4; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}

void cyclic_karatsuba8_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[8], src2x16[8];
    int16x16_t desx16[8];

    for(size_t i = 0; i < 8; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    cyclic_karatsuba8_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 8; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}

void negacyclic_karatsuba8_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[8], src2x16[8];
    int16x16_t desx16[8];

    for(size_t i = 0; i < 8; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    negacyclic_karatsuba8_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 8; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}

void karatsuba_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, size_t len){

    if(len == 1){
        int16x16_t desx16, src1x16, src2x16;
        src1x16 = load_int16x16((int16x16_t*)src1);
        src2x16 = load_int16x16((int16x16_t*)src2);
        desx16 = montmulmod_int16x16(src1x16, src2x16);
        store_int16x16((int16x16_t*)des, desx16);
        return;
    }

    if(len == 2){

        int16x16_t src1x16[2], src2x16[2];
        int16x16_t desx16[3];

        for(size_t i = 0; i < 2; i++){
            src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
            src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
        }

        karatsuba2_Rmod_int16x16(desx16, src1x16, src2x16);

        for(size_t i = 0; i < 3; i++){
            store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
        }

        return;
    }

    if(len == 4){

        int16x16_t src1x16[4], src2x16[4];
        int16x16_t desx16[7];

        for(size_t i = 0; i < 4; i++){
            src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
            src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
        }

        karatsuba4_Rmod_int16x16(desx16, src1x16, src2x16);

        for(size_t i = 0; i < 7; i++){
            store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
        }

        return;

    }

    if(len == 8){

        int16x16_t src1x16[8], src2x16[8];
        int16x16_t desx16[15];

        for(size_t i = 0; i < 8; i++){
            src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
            src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
        }

        karatsuba8_Rmod_int16x16(desx16, src1x16, src2x16);

        for(size_t i = 0; i < 15; i++){
            store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
        }

        return;

    }

    if(len == 16){

        int16x16_t src1x16[16], src2x16[16];
        int16x16_t desx16[31];

        for(size_t i = 0; i < 16; i++){
            src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
            src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
        }

        karatsuba16_Rmod_int16x16(desx16, src1x16, src2x16);

        for(size_t i = 0; i < 31; i++){
            store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
        }

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

    karatsuba_Rmod_int16x16_(des, src1, src2, len / 2);
    karatsuba_Rmod_int16x16_(des + len * 16, src1 + (len / 2) * 16, src2 + (len / 2) * 16, len / 2);
    karatsuba_Rmod_int16x16_(desmid, src1mid, src2mid, len / 2);

    for(size_t i = 0; i < len - 1; i++){
        sub_int16x16_(desmid + i * 16, desmid + i * 16, des + i * 16);
        sub_int16x16_(desmid + i * 16, desmid + i * 16, des + (i + len) * 16);
    }

    for(size_t i = 0; i < len - 1; i++){
        add_int16x16_(des + (i + (len / 2)) * 16, des + (i + (len / 2)) * 16, desmid + i * 16);
    }

}

// notice that we don't need to scale twiddle
void weighted_Karatsuba_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, size_t len, int16_t *twiddle){

    int16_t buff[(2 * len - 1) * 16];

    karatsuba_Rmod_int16x16_(buff, src1, src2, len);

    for(size_t i = len; i < 2 * len - 1; i++){
        mul_int16x16_(buff + i * 16, buff + i * 16, twiddle);
        add_int16x16_(buff + (i - len) * 16, buff + (i - len) * 16, buff + i * 16);
    }

    memmove(des, buff, len * 16 * sizeof(int16_t));

}







