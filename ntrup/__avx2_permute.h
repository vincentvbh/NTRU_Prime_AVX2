#ifndef __AVX2_PERMUTE_H
#define __AVX2_PERMUTE_H

#include <stdint.h>

extern void __asm_permute_test(int16_t *des, int16_t *src);

void transpose(int16_t out[16][16],const int16_t in[16][16]);

void transpose_partial3(int16_t out[16][16], int16_t *in0, int16_t *in1, int16_t *in2, int16_t *in3, int16_t *int4, int16_t *in5);
void transpose_partial3_inv(int16_t *out0, int16_t *out1, int16_t *out2, int16_t *out3, int16_t *out4, int16_t *out5 ,const int16_t in[16][16]);


#endif

