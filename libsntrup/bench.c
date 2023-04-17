
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
#include "avx.h"

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

    for(size_t i = p; i < ARRAY_N; i++){
        poly1[i] = poly2[i] = 0;
    }

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        mult768_over64(res, poly1, poly2);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("polymul cycles: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = cpucycles();
        ntrup_mul(res, poly1, poly2);
        t1 = cpucycles();
        cycles[i] = t1 - t0;
    }
    qsort(cycles, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("ntrup mul cycles: ");
    printf("%lu, %lu, %lu\n",
            cycles[ITERATIONS >> 2], cycles[ITERATIONS >> 1], cycles[(ITERATIONS >> 1) + (ITERATIONS >> 2)]);



}




