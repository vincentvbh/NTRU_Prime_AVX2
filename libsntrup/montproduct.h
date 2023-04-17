#ifndef montproduct_H
#define montproduct_H


#include "stdint.h"
#include "immintrin.h"

// note : montproduct with v4158_16 to multiply by R=2^16

#define q 4591
#define qinv  15631    /* reciprocal of q mod 2^16 */
#define Rx64   2721    /* (2**16)*64 mod q  */
#define q15 7          /* round(2^15/q)     */

static inline
__m256i montproduct( __m256i x , __m256i y )
{
  __m256i lo = _mm256_mullo_epi16( x , y );
  __m256i hi = _mm256_mulhi_epi16( x , y );
  __m256i d = _mm256_mullo_epi16( lo , _mm256_set1_epi16(qinv) );
  __m256i e = _mm256_mulhi_epi16( d , _mm256_set1_epi16(q) );
  return _mm256_sub_epi16(hi,e);
}


static inline
void montproduct_Rx64_p16( int16_t* r, const int16_t* x )
{
  __m256i xx = _mm256_loadu_si256( (__m256i const *) x );
  __m256i rr = montproduct( xx , _mm256_set1_epi16(Rx64) );
  _mm256_storeu_si256( (__m256i *)r , rr );
}




#endif
