
#include "__avx2_basemul_FFT.h"

void cyclic_FFT2_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[2], src2x16[2];
    int16x16_t desx16[2];

    for(size_t i = 0; i < 2; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    cyclic_FFT2_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 2; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}



void cyclic_FFT4_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[4], src2x16[4];
    int16x16_t desx16[4];

    for(size_t i = 0; i < 4; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    cyclic_FFT4_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 4; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}

void cyclic_FFT8_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){

    int16x16_t src1x16[8], src2x16[8];
    int16x16_t desx16[8];

    for(size_t i = 0; i < 8; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    cyclic_FFT8_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 8; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}



void cyclic_FFT16_Rmod_int16x16_(int16_t *des, int16_t *src1, int16_t *src2){


    int16x16_t src1x16[16], src2x16[16];
    int16x16_t desx16[16];

    for(size_t i = 0; i < 16; i++){
        src1x16[i] = load_int16x16((int16x16_t*)(src1 + i * 16));
        src2x16[i] = load_int16x16((int16x16_t*)(src2 + i * 16));
    }

    cyclic_FFT16_Rmod_int16x16(desx16, src1x16, src2x16);

    for(size_t i = 0; i < 16; i++){
        store_int16x16((int16x16_t*)(des + i * 16), desx16[i]);
    }

}



// TODO: don't modify inputs
void cyclic_FFT16_precompute_int16x16_(
    int16_t *des_ptr, int16_t *src1_ptr, int16_t *src2_ptr, int16_t *src2inv_ptr){


    int16x16_t src1[16], src2[16], src2inv[16];
    int16x16_t buff[16];
    int16x16_t t1[16];

    src1[ 0] = load_int16x16((int16x16_t*)(src1_ptr +  0 * 16));
    src1[ 2] = load_int16x16((int16x16_t*)(src1_ptr +  2 * 16));
    src1[ 4] = load_int16x16((int16x16_t*)(src1_ptr +  4 * 16));
    src1[ 6] = load_int16x16((int16x16_t*)(src1_ptr +  6 * 16));
    src1[ 8] = load_int16x16((int16x16_t*)(src1_ptr +  8 * 16));
    src1[10] = load_int16x16((int16x16_t*)(src1_ptr + 10 * 16));
    src1[12] = load_int16x16((int16x16_t*)(src1_ptr + 12 * 16));
    src1[14] = load_int16x16((int16x16_t*)(src1_ptr + 14 * 16));

    t1[ 0] = add_int16x16(src1[ 0], src1[ 8]);
    t1[ 8] = sub_int16x16(src1[ 0], src1[ 8]);
    t1[ 4] = add_int16x16(src1[ 4], src1[12]);
    t1[12] = sub_int16x16(src1[ 4], src1[12]);

    src1[ 0] = barrett_int16x16(add_int16x16(t1[ 0], t1[ 4]));
    src1[ 4] = barrett_int16x16(sub_int16x16(t1[ 0], t1[ 4]));

    for(size_t i = 0; i < 4; i++){
        src2[i] = load_int16x16((int16x16_t*)(src2_ptr + i * 16));
        src2inv[i] = load_int16x16((int16x16_t*)(src2inv_ptr + i * 16));
    }

    t1[ 2] = add_int16x16(src1[ 2], src1[10]);
    t1[10] = sub_int16x16(src1[ 2], src1[10]);
    t1[ 6] = add_int16x16(src1[ 6], src1[14]);
    t1[14] = sub_int16x16(src1[ 6], src1[14]);

    src1[ 2] = barrett_int16x16(add_int16x16(t1[ 2], t1[ 6]));
    src1[ 6] = barrett_int16x16(sub_int16x16(t1[ 2], t1[ 6]));

    t1[ 0] = add_int16x16(src1[ 0], src1[ 2]);
    t1[ 2] = sub_int16x16(src1[ 0], src1[ 2]);

    src1[ 1] = load_int16x16((int16x16_t*)(src1_ptr +  1 * 16));
    src1[ 3] = load_int16x16((int16x16_t*)(src1_ptr +  3 * 16));
    src1[ 5] = load_int16x16((int16x16_t*)(src1_ptr +  5 * 16));
    src1[ 7] = load_int16x16((int16x16_t*)(src1_ptr +  7 * 16));
    src1[ 9] = load_int16x16((int16x16_t*)(src1_ptr +  9 * 16));
    src1[11] = load_int16x16((int16x16_t*)(src1_ptr + 11 * 16));
    src1[13] = load_int16x16((int16x16_t*)(src1_ptr + 13 * 16));
    src1[15] = load_int16x16((int16x16_t*)(src1_ptr + 15 * 16));

    t1[ 1] = add_int16x16(src1[ 1], src1[ 9]);
    t1[ 9] = sub_int16x16(src1[ 1], src1[ 9]);
    t1[ 5] = add_int16x16(src1[ 5], src1[13]);
    t1[13] = sub_int16x16(src1[ 5], src1[13]);

    src1[ 1] = barrett_int16x16(add_int16x16(t1[ 1], t1[ 5]));
    src1[ 5] = barrett_int16x16(sub_int16x16(t1[ 1], t1[ 5]));

    for(size_t i = 4; i < 8; i++){
        src2[i] = load_int16x16((int16x16_t*)(src2_ptr + i * 16));
        src2inv[i] = load_int16x16((int16x16_t*)(src2inv_ptr + i * 16));
    }

    t1[ 3] = add_int16x16(src1[ 3], src1[11]);
    t1[11] = sub_int16x16(src1[ 3], src1[11]);
    t1[ 7] = add_int16x16(src1[ 7], src1[15]);
    t1[15] = sub_int16x16(src1[ 7], src1[15]);

    src1[ 3] = barrett_int16x16(add_int16x16(t1[ 3], t1[ 7]));
    src1[ 7] = barrett_int16x16(sub_int16x16(t1[ 3], t1[ 7]));

    for(size_t i = 8; i < 12; i++){
        src2[i] = load_int16x16((int16x16_t*)(src2_ptr + i * 16));
        src2inv[i] = load_int16x16((int16x16_t*)(src2inv_ptr + i * 16));
    }

    t1[ 1] = add_int16x16(src1[ 1], src1[ 3]);
    t1[ 3] = sub_int16x16(src1[ 1], src1[ 3]);

    src1[0] = add_int16x16(t1[0], t1[1]);
    src1[1] = sub_int16x16(t1[0], t1[1]);

    negacyclic_schoolbook2_precompute_Rmod_int16x16(buff + 2, t1 + 2, src2 + 2, src2inv + 2);

    t1[0] = montmulmod_precompute_int16x16(src1[0], src2[0], src2inv[0]);
    t1[1] = montmulmod_precompute_int16x16(src1[1], src2[1], src2inv[1]);

    buff[0] = add_int16x16(t1[0], t1[1]);
    buff[1] = sub_int16x16(t1[0], t1[1]);

    for(size_t i = 12; i < 16; i++){
        src2[i] = load_int16x16((int16x16_t*)(src2_ptr + i * 16));
        src2inv[i] = load_int16x16((int16x16_t*)(src2inv_ptr + i * 16));
    }

    negacyclic_karatsuba4_precompute_Rmod_int16x16(t1 + 4, src1 + 4, src2 + 4, src2inv + 4);

    t1[0] = barrett_int16x16(add_int16x16(buff[0], buff[2]));
    t1[2] = barrett_int16x16(sub_int16x16(buff[0], buff[2]));
    t1[1] = barrett_int16x16(add_int16x16(buff[1], buff[3]));
    t1[3] = barrett_int16x16(sub_int16x16(buff[1], buff[3]));

    for(size_t i = 0; i < 4; i++){
        buff[i + 4] = sub_int16x16(t1[i], t1[i + 4]);
        buff[i    ] = add_int16x16(t1[i], t1[i + 4]);
    }

    negacyclic_karatsuba8_precompute_Rmod_int16x16(buff + 8, t1 + 8, src2 + 8, src2inv + 8);

    for(size_t i = 0; i < 8; i++){
        store_int16x16((int16x16_t*)(des_ptr +       i * 16), add_int16x16(buff[i], buff[i + 8]));
        store_int16x16((int16x16_t*)(des_ptr + (i + 8) * 16), sub_int16x16(buff[i], buff[i + 8]));
    }

}

void weighted16_Rmod_FFT2_int16x16(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_CT_table, int16_t *_GS_table){

    int16_t buff1[16 * 16], buff2[16 * 16];
    int16_t buff3[32 * 16];
    int16_t k1[8 * 16], k2[8 * 16];
    int16_t k3[48 * 16];

    int16x16_t a0, a1, a2, a3, a4, a5, a6, a7;
    int16x16_t t0, t1;
    int16x16_t twiddle, twiddle_inv;

    (void)a5; (void)a7;

    twiddle = load_int16x16((int16x16_t*)_CT_table);
    twiddle_inv = mullo_int16x16(twiddle, dup_const_int16x16(Qinv));


    for(size_t i = 0; i < 4; i++){
        a0 = load_int16x16((int16x16_t*)(src1 + (i +  0) * 16));
        a1 = load_int16x16((int16x16_t*)(src1 + (i +  4) * 16));
        a2 = load_int16x16((int16x16_t*)(src1 + (i +  8) * 16));
        a3 = load_int16x16((int16x16_t*)(src1 + (i + 12) * 16));
        t0 = montmulmod_precompute_int16x16(a2, twiddle, twiddle_inv);
        t1 = montmulmod_precompute_int16x16(a3, twiddle, twiddle_inv);
        a2 = sub_int16x16(a0, t0);
        a3 = sub_int16x16(a1, t1);
        a0 = add_int16x16(a0, t0);
        a1 = add_int16x16(a1, t1);
        store_int16x16((int16x16_t*)(k1 + (i + 0) * 16), add_int16x16(a0, a1));
        store_int16x16((int16x16_t*)(k1 + (i + 4) * 16), add_int16x16(a2, a3));
        store_int16x16((int16x16_t*)(buff1 + (i +  0) * 16), a0);
        store_int16x16((int16x16_t*)(buff1 + (i +  4) * 16), a1);
        store_int16x16((int16x16_t*)(buff1 + (i +  8) * 16), a2);
        store_int16x16((int16x16_t*)(buff1 + (i + 12) * 16), a3);
    }

    for(size_t i = 0; i < 4; i++){
        a0 = load_int16x16((int16x16_t*)(src2 + (i +  0) * 16));
        a1 = load_int16x16((int16x16_t*)(src2 + (i +  4) * 16));
        a2 = load_int16x16((int16x16_t*)(src2 + (i +  8) * 16));
        a3 = load_int16x16((int16x16_t*)(src2 + (i + 12) * 16));
        t0 = montmulmod_precompute_int16x16(a2, twiddle, twiddle_inv);
        t1 = montmulmod_precompute_int16x16(a3, twiddle, twiddle_inv);
        a2 = sub_int16x16(a0, t0);
        a3 = sub_int16x16(a1, t1);
        a0 = add_int16x16(a0, t0);
        a1 = add_int16x16(a1, t1);
        store_int16x16((int16x16_t*)(k2 + (i + 0) * 16), add_int16x16(a0, a1));
        store_int16x16((int16x16_t*)(k2 + (i + 4) * 16), add_int16x16(a2, a3));
        store_int16x16((int16x16_t*)(buff2 + (i +  0) * 16), a0);
        store_int16x16((int16x16_t*)(buff2 + (i +  4) * 16), a1);
        store_int16x16((int16x16_t*)(buff2 + (i +  8) * 16), a2);
        store_int16x16((int16x16_t*)(buff2 + (i + 12) * 16), a3);
    }

    __asm_schoolbook4(buff3 + 0 * 16, buff1 + 0 * 16, buff2 + 0 * 16, const_buff);
    __asm_schoolbook4(buff3 + 8 * 16, buff1 + 4 * 16, buff2 + 4 * 16, const_buff);
    __asm_schoolbook4(k3 + 0 * 16, k1 + 0 * 16, k2 + 0 * 16, const_buff);

    store_int16x16((int16x16_t*)(buff3 + 7 * 16), barrett_int16x16(sub_int16x16(sub_int16x16(
                                                        load_int16x16((int16x16_t*)(k3 + 3 * 16)),
                                                            load_int16x16((int16x16_t*)(buff3 + 3 * 16))),
                                                            load_int16x16((int16x16_t*)(buff3 + 11 * 16)))));

    for(size_t i = 0; i < 3; i++){

        a2 = barrett_int16x16(sub_int16x16(sub_int16x16(
                                load_int16x16((int16x16_t*)(k3 + (i + 0) * 16)),
                                    load_int16x16((int16x16_t*)(buff3 + (i + 0) * 16))),
                                    load_int16x16((int16x16_t*)(buff3 + (i + 8) * 16))));
        a4 = barrett_int16x16(sub_int16x16(sub_int16x16(
                                load_int16x16((int16x16_t*)(k3 + (i + 4) * 16)),
                                    load_int16x16((int16x16_t*)(buff3 + (i + 4) * 16))),
                                    load_int16x16((int16x16_t*)(buff3 + (i + 12) * 16))));
        a2 = add_int16x16(load_int16x16((int16x16_t*)(buff3 + (i + 4) * 16)), a2);
        a4 = add_int16x16(load_int16x16((int16x16_t*)(buff3 + (i + 8) * 16)), a4);

        a0 = load_int16x16((int16x16_t*)(buff3 + (i + 0) * 16));
        a6 = load_int16x16((int16x16_t*)(buff3 + (i + 12) * 16));

        a4 = montmulmod_precompute_int16x16(a4, twiddle, twiddle_inv);
        a6 = montmulmod_precompute_int16x16(a6, twiddle, twiddle_inv);

        a0 = add_int16x16(a0, a4);
        a2 = add_int16x16(a2, a6);

        store_int16x16((int16x16_t*)(buff3 + (i + 0) * 16), a0);
        store_int16x16((int16x16_t*)(buff3 + (i + 4) * 16), a2);

    }

    a0 = load_int16x16((int16x16_t*)(buff3 +  3 * 16));
    a2 = load_int16x16((int16x16_t*)(buff3 + 11 * 16));
    t0 = montmulmod_precompute_int16x16(a2, twiddle, twiddle_inv);
    a0 = add_int16x16(a0, t0);
    store_int16x16((int16x16_t*)(buff3 + 3 * 16), a0);

    __asm_schoolbook4(buff3 + 8 * 16, buff1 + 8 * 16, buff2 + 8 * 16, const_buff);
    __asm_schoolbook4(buff3 + 16 * 16, buff1 + 12 * 16, buff2 + 12 * 16, const_buff);
    __asm_schoolbook4(k3 + 8 * 16, k1 + 4 * 16, k2 + 4 * 16, const_buff);

    store_int16x16((int16x16_t*)(buff3 + 15 * 16), barrett_int16x16(sub_int16x16(sub_int16x16(
                                                        load_int16x16((int16x16_t*)(k3 + 11 * 16)),
                                                            load_int16x16((int16x16_t*)(buff3 + 11 * 16))),
                                                            load_int16x16((int16x16_t*)(buff3 + 19 * 16)))));

    for(size_t i = 0; i < 3; i++){

        a2 = barrett_int16x16(sub_int16x16(sub_int16x16(
                                load_int16x16((int16x16_t*)(k3 + (i + 8) * 16)),
                                    load_int16x16((int16x16_t*)(buff3 + (i + 8) * 16))),
                                    load_int16x16((int16x16_t*)(buff3 + (i + 16) * 16))));
        a4 = barrett_int16x16(sub_int16x16(sub_int16x16(
                                load_int16x16((int16x16_t*)(k3 + (i + 12) * 16)),
                                    load_int16x16((int16x16_t*)(buff3 + (i + 12) * 16))),
                                    load_int16x16((int16x16_t*)(buff3 + (i + 20) * 16))));
        a2 = add_int16x16(load_int16x16((int16x16_t*)(buff3 + (i + 12) * 16)), a2);
        a4 = add_int16x16(load_int16x16((int16x16_t*)(buff3 + (i + 16) * 16)), a4);

        a0 = load_int16x16((int16x16_t*)(buff3 + (i + 8) * 16));
        a6 = load_int16x16((int16x16_t*)(buff3 + (i + 20) * 16));

        a4 = montmulmod_precompute_int16x16(a4, twiddle, twiddle_inv);
        a6 = montmulmod_precompute_int16x16(a6, twiddle, twiddle_inv);

        a0 = sub_int16x16(a0, a4);
        a2 = sub_int16x16(a2, a6);

        store_int16x16((int16x16_t*)(buff3 + (i + 8) * 16), a0);
        store_int16x16((int16x16_t*)(buff3 + (i + 12) * 16), a2);

    }

    a0 = load_int16x16((int16x16_t*)(buff3 + 11 * 16));
    a2 = load_int16x16((int16x16_t*)(buff3 + 19 * 16));
    t0 = montmulmod_precompute_int16x16(a2, twiddle, twiddle_inv);
    a0 = sub_int16x16(a0, t0);
    store_int16x16((int16x16_t*)(buff3 + 11 * 16), a0);

    twiddle = load_int16x16((int16x16_t*)_GS_table);
    twiddle_inv = mullo_int16x16(twiddle, dup_const_int16x16(Qinv));


    for(size_t i = 0; i < 8; i++){
        a0 = load_int16x16((int16x16_t*)(buff3 + (i + 0) * 16));
        a2 = load_int16x16((int16x16_t*)(buff3 + (i + 8) * 16));
        t0 = sub_int16x16(a0, a2);
        a0 = barrett_int16x16(add_int16x16(a0, a2));
        a2 = montmulmod_precompute_int16x16(t0, twiddle, twiddle_inv);
        store_int16x16((int16x16_t*)(des + (i + 0) * 16), a0);
        store_int16x16((int16x16_t*)(des + (i + 8) * 16), a2);
    }

}

