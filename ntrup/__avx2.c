
#include "__avx2.h"

void polymul(int16_t *des, const int16_t *src1, const int16_t *src2){

    int16_t src1_NTT[ARRAY_N];
    int16_t src2_NTT[ARRAY_N];
    int16_t des_NTT[(ARRAY_N / 2)];

    int16_t buff1[256], buff2[256], buff3[256];

    // for(size_t i = 0; i < ARRAY_N; i += 16){
    //     t = i / 16;
    //     memmove(src1_NTT + ( (t % 2) * 51 +  (t % 17) * 3 + (t % 3)) * 16, src1 + i, 16 * sizeof(int16_t));
    // }

    for(size_t i = 0; i < ARRAY_N; i += (ARRAY_N / 2)){
        for(size_t j = 0; j < (ARRAY_N / (2 * 17)); j += 16){
            __asm_rader17_pre_tbl(src1_NTT + i + j, src1, const_twiddle17, const_twiddle17pre, const_buff,
                                  permute_tbl + (i / (ARRAY_N / 2)) * 51 + (j / 16) * 17);
            __asm_rader17_pre_tbl(src2_NTT + i + j, src2, const_twiddle17, const_twiddle17pre, const_buff,
                                  permute_tbl + (i / (ARRAY_N / 2)) * 51 + (j / 16) * 17);
        }
    }

    __asm_3x2_pre(src1_NTT, src1_NTT + 17 * 48, OMEGA3_buff, const_buff);
    __asm_3x2_pre(src2_NTT, src2_NTT + 17 * 48, OMEGA3_buff, const_buff);

// ================================
// Cooley--Tukey

    for(size_t i = 0; i < 48; i += 16){
        transpose((int16_t (*)[16])(src1_NTT + i * 16),
                  (int16_t (*)[16])(src1_NTT + i * 16));
        transpose((int16_t (*)[16])(src2_NTT + i * 16),
                  (int16_t (*)[16])(src2_NTT + i * 16));
    }

    weighted16_Rmod_FFT2_int16x16(src1_NTT, src1_NTT, src2_NTT, twiddle_CT, twiddle_GS);
    weighted16_Rmod_FFT2_int16x16(src1_NTT + 16 * 16, src1_NTT + 16 * 16, src2_NTT + 16 * 16, twiddle_CT + 16, twiddle_GS + 16);
    weighted16_Rmod_FFT2_int16x16(src1_NTT + 32 * 16, src1_NTT + 32 * 16, src2_NTT + 32 * 16, twiddle_CT + 32, twiddle_GS + 32);

    for(size_t i = 0; i < 48; i += 16){
        transpose((int16_t (*)[16])(src1_NTT + i * 16),
                  (int16_t (*)[16])(src1_NTT + i * 16));
    }

// ================================
// Karatsuba

    for(size_t i = 0; i < 48; i += 16){
        transpose((int16_t (*)[16])(src1_NTT + (ARRAY_N / 2) + i * 16),
                  (int16_t (*)[16])(src1_NTT + (ARRAY_N / 2) + i * 16));
        transpose((int16_t (*)[16])(src2_NTT + (ARRAY_N / 2) + i * 16),
                  (int16_t (*)[16])(src2_NTT + (ARRAY_N / 2) + i * 16));
    }

    for(size_t i = 0; i < 48; i += 16){
        __asm_weighted_double_karatsuba16(des_NTT + (ARRAY_N / 2) + i * 16 - (ARRAY_N / 2),
                           src1_NTT + (ARRAY_N / 2) + i * 16, src2_NTT + (ARRAY_N / 2) + i * 16,
                           const_buff,
                           const_twiddle51n + i);
    }

    for(size_t i = 0; i < 48; i += 16){
        transpose((int16_t (*)[16])(src1_NTT + (ARRAY_N / 2) + i * 16),
                  (int16_t (*)[16])(des_NTT + (ARRAY_N / 2) + i * 16 - (ARRAY_N / 2)));
    }

// ================================
// remainings

    transpose_partial3((int16_t (*)[16])buff1, src1_NTT                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                               src1_NTT + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                               src1_NTT                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                               src1_NTT + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                               src1_NTT                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                               src1_NTT + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));
    transpose_partial3((int16_t (*)[16])buff2, src2_NTT                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                                               src2_NTT + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                                               src2_NTT                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                                               src2_NTT + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                                               src2_NTT                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                                               src2_NTT + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)));

    __asm_weighted_double_karatsuba16(buff3, buff1, buff2, const_buff, const_twiddle102_last);

    transpose_partial3_inv(src1_NTT                 + 48 * (ARRAY_N / (2 * 17 * 3)),
                           src1_NTT + (ARRAY_N / 2) + 48 * (ARRAY_N / (2 * 17 * 3)),
                           src1_NTT                 + 49 * (ARRAY_N / (2 * 17 * 3)),
                           src1_NTT + (ARRAY_N / 2) + 49 * (ARRAY_N / (2 * 17 * 3)),
                           src1_NTT                 + 50 * (ARRAY_N / (2 * 17 * 3)),
                           src1_NTT + (ARRAY_N / 2) + 50 * (ARRAY_N / (2 * 17 * 3)), (int16_t (*)[16])buff3);

    __asm_3x2_post(src1_NTT, src1_NTT + 17 * 48, OMEGA3INV_buff, const_buff);

    for(size_t i = 0; i < ARRAY_N; i += (ARRAY_N / 2)){
        for(size_t j = 0; j < (ARRAY_N / (2 * 17)); j += 16){
            __asm_rader17_post_tbl(des, src1_NTT + i + j, const_twiddle17inv, const_twiddle17invpre, const_buff,
                                   permute_tbl + (i / (ARRAY_N / 2)) * 51 + (j / 16) * 17);
        }
    }


}

void ntrup_mul(int16_t *des, const int16_t *src1, const int16_t *src2){

    int16_t buff[1632];
    int16x16_t t;

    size_t i;

    polymul(buff, src1, src2);

    for(i = 0; i + 16 < p - 1; i += 16){
        store_int16x16((int16x16_t*)(buff + i), add_int16x16(load_int16x16((int16x16_t*)(buff + i)),
                                                             load_int16x16((int16x16_t*)(buff + i + p))));
    }
    for(; i < p - 1; i++){
        buff[i] = buff[i] + buff[i + p];
    }

    for(i = 1; i + 16 < p; i += 16){
        store_int16x16((int16x16_t*)(buff + i), add_int16x16(load_int16x16((int16x16_t*)(buff + i)),
                                                             load_int16x16((int16x16_t*)(buff + i + p - 1))));
    }

    for(; i < p; i++){
        buff[i] = buff[i] + buff[i + p - 1];
    }

    t = dup_const_int16x16(FINAL_SCALE_Rmod);

    for(i = 0; i < p; i+= 16) {
        store_int16x16((int16x16_t*)(des + i), montmulmod_int16x16(load_int16x16((int16x16_t*)(buff + i)), t));
    }

}

