#ifndef avx_H
#define avx_H


#if defined(CRYPTO_NAMESPACE)
#define mult768_over64   CRYPTO_NAMESPACE(mult768_over64)
#endif

#include <inttypes.h>
#include <immintrin.h>

typedef int16_t int16;
typedef int32_t int32;

typedef __m128i int16x8;
typedef __m256i int16x16;

void mult768_over64(int16 *,const int16 *,const int16 *);
void ntrup_mul(int16_t *des, const int16_t *src1, const int16_t *src2);


#endif
