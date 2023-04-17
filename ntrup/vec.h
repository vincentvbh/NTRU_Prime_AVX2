#ifndef VEC_H
#define VEC_H

#include <stdint.h>
#include <stdint.h>
#include <stddef.h>

#include "params.h"


void add_int16x16_(int16_t des[16], int16_t src1[16], int16_t src2[16]);
void sub_int16x16_(int16_t des[16], int16_t src1[16], int16_t src2[16]);
void mul_int16x16_(int16_t des[16], int16_t src1[16], int16_t src2[16]);

void karatsuba_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, size_t len);
void weighted_int16x16_(int16_t *des, int16_t *src1, int16_t *src2, size_t len, int16_t *twiddle);

#endif

