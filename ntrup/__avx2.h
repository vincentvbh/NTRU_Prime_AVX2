#ifndef __AVX2_H
#define __AVX2_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "params.h"
#include "__avx2_wrap.h"
#include "__avx2_permute.h"
#include "__avx2_basemul.h"
#include "__avx2_basemul_Karatsuba.h"
#include "__avx2_basemul_FFT.h"
#include "__avx2_FFT.h"

void polymul(int16_t *des, const int16_t *src1, const int16_t *src2);
void ntrup_mul(int16_t *des, const int16_t *src1, const int16_t *src2);

#endif

