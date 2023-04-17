
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
#include "rader.h"


int16_t id_Q[Q], sqrt_Q[Q];

static int16_t get_sqrt(int16_t i){

    int16_t t;

    coeff_ring.memberZ(&t, &i);

    if(t < 0){
        t += Q;
    }

    t = sqrt_Q[t];

    coeff_ring.memberZ(&t, &t);

    return t;

}

int main(void){

    int16_t twiddle17[17];
    int16_t twiddle3[3];
    int16_t twiddle51p[51];
    int16_t twiddle51n[51];

    int16_t twiddle17_Rmod[17];
    int16_t twiddle17_permuted_Rmod[17];
    int16_t twiddle17_precompute[17];
    int16_t twiddle17_precompute_inv[17];

    int16_t twiddlebuff[17];

    int16_t scale, omega, t;

    for(size_t i = 0; i < Q; i++){
        id_Q[i] = i;
    }

    for(size_t i = 0; i < Q; i++){
        sqrt_Q[i] = -1;
    }

    for(size_t i = 0; i < Q; i++){
        coeff_ring.expZ(&t, id_Q + i, 2);
        if(t < 0){
            t += Q;
        }
        if(sqrt_Q[t] == -1){
            sqrt_Q[t] = i;
        }
    }

    uint64_t tbl_buff[2 * 17 * 3];

    for(size_t i = 0; i < ARRAY_N; i += 16){
        t = i / 16;
        tbl_buff[(t % 2) * 51 +  (t % 17) * 3 + (t % 3)] = i * 2;
    }

    for(size_t i = 0; i < 2; i++){
        for(size_t j = 0; j < 3; j++){
            for(size_t k = 0; k < 17; k++){
                permute_tbl[i * 51 + j * 17 + k] = tbl_buff[i * 51 + j + k * 3];
            }
        }
    }

    omega = OMEGA17;
    getExpVec(twiddle17, omega, 17);

    scale = RmodQ;
    for(size_t i = 0; i < 17; i++){
        coeff_ring.mulZ(twiddle17_Rmod + i, twiddle17 + i, &scale);
    }

    for(size_t i = 1; i < 17; i++){
        twiddle17_permuted_Rmod[i] = twiddle17_Rmod[rader_in_permute[17 - i]];
    }

    for(size_t i = 1; i < 17; i++){
        twiddle17_precompute[i] = twiddle17_permuted_Rmod[i];
    }

    for(size_t i = 1; i < 9; i++){
        coeff_ring.addZ(twiddlebuff + i + 0, twiddle17_precompute + i + 0, twiddle17_precompute + i + 8);
        coeff_ring.subZ(twiddlebuff + i + 8, twiddle17_precompute + i + 0, twiddle17_precompute + i + 8);
    }
    for(size_t i = 1; i < 17; i++){
        twiddle17_precompute[i] = twiddlebuff[i];
    }

    for(size_t i = 1; i < 5; i++){
        coeff_ring.addZ(twiddlebuff + i + 0, twiddle17_precompute + i + 0, twiddle17_precompute + i + 4);
        coeff_ring.subZ(twiddlebuff + i + 4, twiddle17_precompute + i + 0, twiddle17_precompute + i + 4);
    }
    for(size_t i = 1; i < 9; i++){
        twiddle17_precompute[i] = twiddlebuff[i];
    }

    for(size_t i = 1; i < 3; i++){
        coeff_ring.addZ(twiddlebuff + i + 0, twiddle17_precompute + i + 0, twiddle17_precompute + i + 2);
        coeff_ring.subZ(twiddlebuff + i + 2, twiddle17_precompute + i + 0, twiddle17_precompute + i + 2);
    }
    for(size_t i = 1; i < 5; i++){
        twiddle17_precompute[i] = twiddlebuff[i];
    }

    for(size_t i = 1; i < 2; i++){
        coeff_ring.addZ(twiddlebuff + i + 0, twiddle17_precompute + i + 0, twiddle17_precompute + i + 1);
        coeff_ring.subZ(twiddlebuff + i + 1, twiddle17_precompute + i + 0, twiddle17_precompute + i + 1);
    }
    for(size_t i = 1; i < 3; i++){
        twiddle17_precompute[i] = twiddlebuff[i];
    }

    scale = 2;
    for(size_t i = 3; i < 5; i++){
        coeff_ring.mulZ(twiddle17_precompute + i, twiddle17_precompute + i, &scale);
    }

    scale = 4;
    for(size_t i = 5; i < 9; i++){
        coeff_ring.mulZ(twiddle17_precompute + i, twiddle17_precompute + i, &scale);
    }

    scale = 8;
    for(size_t i = 9; i < 17; i++){
        coeff_ring.mulZ(twiddle17_precompute + i, twiddle17_precompute + i, &scale);
    }

    scale = INV16;
    for(size_t i = 1; i < 17; i++){
        coeff_ring.mulZ(twiddle17_precompute + i, twiddle17_precompute + i, &scale);
    }

    for(size_t i = 1; i < 17; i++){
        twiddle17_precompute_inv[i] = twiddle17_precompute[i] * Qinv;
    }

    for(size_t i = 1; i < 17; i++){
        for(size_t j = 0; j < 16; j++){
            assert(const_twiddle17[i * 16 + j] == twiddle17_precompute[i]);
            assert(const_twiddle17pre[i * 16 + j] == twiddle17_precompute_inv[i]);
        }
    }

    omega = OMEGA3;
    getExpVec(twiddle3, omega, 3);

    for(size_t i = 0; i < 17; i++){
        for(size_t j = 0; j < 3; j++){
            coeff_ring.mulZ(twiddle51p + i * 3 + j, twiddle17 + i, twiddle3 + j);
        }
    }
    for(size_t i = 0; i < 51; i++){
        twiddle51n[i] = -twiddle51p[i];
    }

    for(size_t i = 0; i < 48; i++){
        twiddle_CT[i] = get_sqrt(twiddle51p[i]);
        coeff_ring.expZ(twiddle_GS + i, twiddle_CT + i, Q - 2);
    }

    scale = RmodQ;
    for(size_t i = 0; i < 48; i++){
        coeff_ring.mulZ(twiddle_CT + i, twiddle_CT + i, &scale);
        coeff_ring.mulZ(twiddle_GS + i, twiddle_GS + i, &scale);
    }

    twiddlebuff[0] = twiddle51p[48];
    twiddlebuff[1] = twiddle51n[48];
    twiddlebuff[2] = twiddle51p[49];
    twiddlebuff[3] = twiddle51n[49];
    twiddlebuff[4] = twiddle51p[50];
    twiddlebuff[5] = twiddle51n[50];

    scale = RmodQ;
    for(size_t i = 0; i < 51; i++){
        coeff_ring.mulZ(twiddle51p + i, twiddle51p + i, &scale);
        coeff_ring.mulZ(twiddle51n + i, twiddle51n + i, &scale);
    }
    for(size_t i = 0; i < 6; i++){
        coeff_ring.mulZ(twiddlebuff + i, twiddlebuff + i, &scale);
    }

    for(size_t i = 0; i < 48; i++){
        assert(twiddle51p[i] == const_twiddle51p[i]);
    }
    for(size_t i = 0; i < 48; i++){
        assert(twiddle51n[i] == const_twiddle51n[i]);
    }
    for(size_t i = 0; i < 6; i++){
        assert(twiddlebuff[i] == const_twiddle102_last[i]);
    }

    omega = OMEGA17INV;
    getExpVec(twiddle17, omega, 17);

    scale = RmodQ;
    for(size_t i = 0; i < 17; i++){
        coeff_ring.mulZ(twiddle17_Rmod + i, twiddle17 + i, &scale);
    }

    for(size_t i = 1; i < 17; i++){
        twiddle17_permuted_Rmod[i] = twiddle17_Rmod[rader_in_permute[17 - i]];
    }

    for(size_t i = 1; i < 17; i++){
        twiddle17_precompute[i] = twiddle17_permuted_Rmod[i];
    }

    for(size_t i = 1; i < 9; i++){
        coeff_ring.addZ(twiddlebuff + i + 0, twiddle17_precompute + i + 0, twiddle17_precompute + i + 8);
        coeff_ring.subZ(twiddlebuff + i + 8, twiddle17_precompute + i + 0, twiddle17_precompute + i + 8);
    }
    for(size_t i = 1; i < 17; i++){
        twiddle17_precompute[i] = twiddlebuff[i];
    }

    for(size_t i = 1; i < 5; i++){
        coeff_ring.addZ(twiddlebuff + i + 0, twiddle17_precompute + i + 0, twiddle17_precompute + i + 4);
        coeff_ring.subZ(twiddlebuff + i + 4, twiddle17_precompute + i + 0, twiddle17_precompute + i + 4);
    }
    for(size_t i = 1; i < 9; i++){
        twiddle17_precompute[i] = twiddlebuff[i];
    }

    for(size_t i = 1; i < 3; i++){
        coeff_ring.addZ(twiddlebuff + i + 0, twiddle17_precompute + i + 0, twiddle17_precompute + i + 2);
        coeff_ring.subZ(twiddlebuff + i + 2, twiddle17_precompute + i + 0, twiddle17_precompute + i + 2);
    }
    for(size_t i = 1; i < 5; i++){
        twiddle17_precompute[i] = twiddlebuff[i];
    }

    for(size_t i = 1; i < 2; i++){
        coeff_ring.addZ(twiddlebuff + i + 0, twiddle17_precompute + i + 0, twiddle17_precompute + i + 1);
        coeff_ring.subZ(twiddlebuff + i + 1, twiddle17_precompute + i + 0, twiddle17_precompute + i + 1);
    }
    for(size_t i = 1; i < 3; i++){
        twiddle17_precompute[i] = twiddlebuff[i];
    }

    scale = 2;
    for(size_t i = 3; i < 5; i++){
        coeff_ring.mulZ(twiddle17_precompute + i, twiddle17_precompute + i, &scale);
    }

    scale = 4;
    for(size_t i = 5; i < 9; i++){
        coeff_ring.mulZ(twiddle17_precompute + i, twiddle17_precompute + i, &scale);
    }

    scale = 8;
    for(size_t i = 9; i < 17; i++){
        coeff_ring.mulZ(twiddle17_precompute + i, twiddle17_precompute + i, &scale);
    }

    scale = INV16;
    for(size_t i = 1; i < 17; i++){
        coeff_ring.mulZ(twiddle17_precompute + i, twiddle17_precompute + i, &scale);
    }

    for(size_t i = 1; i < 17; i++){
        twiddle17_precompute_inv[i] = twiddle17_precompute[i] * Qinv;
    }

    for(size_t i = 1; i < 17; i++){
        for(size_t j = 0; j < 16; j++){
            assert(const_twiddle17inv[i * 16 + j] == twiddle17_precompute[i]);
            assert(const_twiddle17invpre[i * 16 + j] == twiddle17_precompute_inv[i]);
        }
    }



    printf("gen finished!\n");



}




