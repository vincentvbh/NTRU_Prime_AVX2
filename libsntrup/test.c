
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

#include "ring.h"
#include "avx.h"

int main(void){

    int16_t poly1[ARRAY_N], poly2[ARRAY_N];
    int16_t ref[ARRAY_N], res[ARRAY_N];

    int16_t scale, twiddle, t;

    for(size_t i = 0; i < p; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

    for(size_t i = p; i < ARRAY_N; i++){
        poly1[i] = poly2[i] = 0;
    }

    for(size_t i = 0; i < p; i++){
        t = rand();
        coeff_ring.memberZ(poly1 + i, &t);
        t = rand();
        coeff_ring.memberZ(poly2 + i, &t);
    }

    for(size_t i = p; i < ARRAY_N; i++){
        poly1[i] = poly2[i] = 0;
    }

    twiddle = 1;
    naive_mulR(ref, poly1, poly2, ARRAY_N, &twiddle, coeff_ring);

    mult768_over64(res, poly1, poly2);

    scale =  FINAL_SCALE;
    for(size_t i = 0; i < ARRAY_N; i++){
        coeff_ring.mulZ(res + i, res + i, &scale);
    }

    for(size_t i = 0; i < ARRAY_N; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("mulcore finished!\n");

    for(size_t i = p + p - 2; i >=p; i--){
        coeff_ring.addZ(ref + i - p, ref + i - p, ref + i);
        coeff_ring.addZ(ref + i - p + 1, ref + i - p + 1, ref + i);
    }

    polymul(res, poly1, poly2);

    for(size_t i = 0; i < p; i++){
        coeff_ring.memberZ(res + i, res + i);
    }

    for(size_t i = 0; i < p; i++){
        if(ref[i] != res[i]){
            printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
        }
    }

    printf("polymul finished!\n");





}




