
#include "__avx2_wrap.h"
#include "__avx2_permute.h"

void transpose(int16_t out[16][16], const int16_t in[16][16])
{
    int16x16_t a0 = load_int16x16((int16x16_t*)&in[0][0]);
    int16x16_t a1 = load_int16x16((int16x16_t*)&in[1][0]);
    int16x16_t b0 = _mm256_unpacklo_epi16(a0,a1);
    int16x16_t b1 = _mm256_unpackhi_epi16(a0,a1);
    int16x16_t a2 = load_int16x16((int16x16_t*)&in[2][0]);
    int16x16_t a3 = load_int16x16((int16x16_t*)&in[3][0]);
    int16x16_t b2 = _mm256_unpacklo_epi16(a2,a3);
    int16x16_t b3 = _mm256_unpackhi_epi16(a2,a3);
    int16x16_t c0 = _mm256_unpacklo_epi32(b0,b2);
    int16x16_t c2 = _mm256_unpackhi_epi32(b0,b2);
    int16x16_t c1 = _mm256_unpacklo_epi32(b1,b3);
    int16x16_t c3 = _mm256_unpackhi_epi32(b1,b3);
    int16x16_t a4 = load_int16x16((int16x16_t*)&in[4][0]);
    int16x16_t a5 = load_int16x16((int16x16_t*)&in[5][0]);
    int16x16_t b4 = _mm256_unpacklo_epi16(a4,a5);
    int16x16_t b5 = _mm256_unpackhi_epi16(a4,a5);
    int16x16_t a6 = load_int16x16((int16x16_t*)&in[6][0]);
    int16x16_t a7 = load_int16x16((int16x16_t*)&in[7][0]);
    int16x16_t b6 = _mm256_unpacklo_epi16(a6,a7);
    int16x16_t b7 = _mm256_unpackhi_epi16(a6,a7);
    int16x16_t c4 = _mm256_unpacklo_epi32(b4,b6);
    int16x16_t c6 = _mm256_unpackhi_epi32(b4,b6);
    int16x16_t c5 = _mm256_unpacklo_epi32(b5,b7);
    int16x16_t c7 = _mm256_unpackhi_epi32(b5,b7);
    int16x16_t a8 = load_int16x16((int16x16_t*)&in[8][0]);
    int16x16_t a9 = load_int16x16((int16x16_t*)&in[9][0]);
    int16x16_t b8 = _mm256_unpacklo_epi16(a8,a9);
    int16x16_t b9 = _mm256_unpackhi_epi16(a8,a9);
    int16x16_t a10 = load_int16x16((int16x16_t*)&in[10][0]);
    int16x16_t a11 = load_int16x16((int16x16_t*)&in[11][0]);
    int16x16_t b10 = _mm256_unpacklo_epi16(a10,a11);
    int16x16_t b11 = _mm256_unpackhi_epi16(a10,a11);
    int16x16_t c8 = _mm256_unpacklo_epi32(b8,b10);
    int16x16_t c10 = _mm256_unpackhi_epi32(b8,b10);
    int16x16_t c9 = _mm256_unpacklo_epi32(b9,b11);
    int16x16_t c11 = _mm256_unpackhi_epi32(b9,b11);
    int16x16_t a12 = load_int16x16((int16x16_t*)&in[12][0]);
    int16x16_t a13 = load_int16x16((int16x16_t*)&in[13][0]);
    int16x16_t b12 = _mm256_unpacklo_epi16(a12,a13);
    int16x16_t b13 = _mm256_unpackhi_epi16(a12,a13);
    int16x16_t a14 = load_int16x16((int16x16_t*)&in[14][0]);
    int16x16_t a15 = load_int16x16((int16x16_t*)&in[15][0]);
    int16x16_t b14 = _mm256_unpacklo_epi16(a14,a15);
    int16x16_t b15 = _mm256_unpackhi_epi16(a14,a15);
    int16x16_t c12 = _mm256_unpacklo_epi32(b12,b14);
    int16x16_t c14 = _mm256_unpackhi_epi32(b12,b14);
    int16x16_t c13 = _mm256_unpacklo_epi32(b13,b15);
    int16x16_t c15 = _mm256_unpackhi_epi32(b13,b15);

    int16x16_t d0 = _mm256_unpacklo_epi64(c0,c4);
    int16x16_t d4 = _mm256_unpackhi_epi64(c0,c4);
    int16x16_t d8 = _mm256_unpacklo_epi64(c8,c12);
    int16x16_t d12 = _mm256_unpackhi_epi64(c8,c12);
    int16x16_t e0 = _mm256_permute2x128_si256(d0,d8,0x20);
    int16x16_t e8 = _mm256_permute2x128_si256(d0,d8,0x31);
    int16x16_t e4 = _mm256_permute2x128_si256(d4,d12,0x20);
    int16x16_t e12 = _mm256_permute2x128_si256(d4,d12,0x31);
    store_int16x16((int16x16_t*)(out[0]),  e0);
    store_int16x16((int16x16_t*)(out[8]),  e8);
    store_int16x16((int16x16_t*)(out[1]),  e4);
    store_int16x16((int16x16_t*)(out[9]), e12);

    int16x16_t d1 = _mm256_unpacklo_epi64(c1,c5);
    int16x16_t d5 = _mm256_unpackhi_epi64(c1,c5);
    int16x16_t d9 = _mm256_unpacklo_epi64(c9,c13);
    int16x16_t d13 = _mm256_unpackhi_epi64(c9,c13);
    int16x16_t e1 = _mm256_permute2x128_si256(d1,d9,0x20);
    int16x16_t e9 = _mm256_permute2x128_si256(d1,d9,0x31);
    int16x16_t e5 = _mm256_permute2x128_si256(d5,d13,0x20);
    int16x16_t e13 = _mm256_permute2x128_si256(d5,d13,0x31);
    store_int16x16((int16x16_t*)(out[ 4]),  e1);
    store_int16x16((int16x16_t*)(out[12]),  e9);
    store_int16x16((int16x16_t*)(out[ 5]),  e5);
    store_int16x16((int16x16_t*)(out[13]), e13);

    int16x16_t d2 = _mm256_unpacklo_epi64(c2,c6);
    int16x16_t d6 = _mm256_unpackhi_epi64(c2,c6);
    int16x16_t d10 = _mm256_unpacklo_epi64(c10,c14);
    int16x16_t d14 = _mm256_unpackhi_epi64(c10,c14);
    int16x16_t e2 = _mm256_permute2x128_si256(d2,d10,0x20);
    int16x16_t e10 = _mm256_permute2x128_si256(d2,d10,0x31);
    int16x16_t e6 = _mm256_permute2x128_si256(d6,d14,0x20);
    int16x16_t e14 = _mm256_permute2x128_si256(d6,d14,0x31);
    store_int16x16((int16x16_t*)(out[ 2]),  e2);
    store_int16x16((int16x16_t*)(out[ 3]),  e6);
    store_int16x16((int16x16_t*)(out[10]), e10);
    store_int16x16((int16x16_t*)(out[11]), e14);

    int16x16_t d3 = _mm256_unpacklo_epi64(c3,c7);
    int16x16_t d7 = _mm256_unpackhi_epi64(c3,c7);
    int16x16_t d11 = _mm256_unpacklo_epi64(c11,c15);
    int16x16_t d15 = _mm256_unpackhi_epi64(c11,c15);
    int16x16_t e3 = _mm256_permute2x128_si256(d3,d11,0x20);
    int16x16_t e11 = _mm256_permute2x128_si256(d3,d11,0x31);
    int16x16_t e7 = _mm256_permute2x128_si256(d7,d15,0x20);
    int16x16_t e15 = _mm256_permute2x128_si256(d7,d15,0x31);
    store_int16x16((int16x16_t*)(out[ 6]),  e3);
    store_int16x16((int16x16_t*)(out[ 7]),  e7);
    store_int16x16((int16x16_t*)(out[14]), e11);
    store_int16x16((int16x16_t*)(out[15]), e15);

}

void transpose_partial3(int16_t out[16][16], int16_t *in0, int16_t *in1, int16_t *in2,
                                             int16_t *in3, int16_t *in4, int16_t *in5)
{

    int16x16_t _zero;

    _zero = dup_const_int16x16(0);

    int16x16_t a0 = load_int16x16((int16x16_t*)in0);
    int16x16_t a1 = load_int16x16((int16x16_t*)in1);
    int16x16_t b0 = _mm256_unpacklo_epi16(a0,a1);
    int16x16_t b1 = _mm256_unpackhi_epi16(a0,a1);
    int16x16_t a2 = load_int16x16((int16x16_t*)in2);
    int16x16_t a3 = load_int16x16((int16x16_t*)in3);
    int16x16_t b2 = _mm256_unpacklo_epi16(a2,a3);
    int16x16_t b3 = _mm256_unpackhi_epi16(a2,a3);
    int16x16_t c0 = _mm256_unpacklo_epi32(b0,b2);
    int16x16_t c2 = _mm256_unpackhi_epi32(b0,b2);
    int16x16_t c1 = _mm256_unpacklo_epi32(b1,b3);
    int16x16_t c3 = _mm256_unpackhi_epi32(b1,b3);
    int16x16_t a4 = load_int16x16((int16x16_t*)in4);
    int16x16_t a5 = load_int16x16((int16x16_t*)in5);
    int16x16_t b4 = _mm256_unpacklo_epi16(a4,a5);
    int16x16_t b5 = _mm256_unpackhi_epi16(a4,a5);
    // int16x16_t a6 = _zero;
    // int16x16_t a7 = _zero;
    // int16x16_t b6 = _zero;
    // int16x16_t b7 = _zero;
    int16x16_t c4 = _mm256_unpacklo_epi32(b4,_zero);
    int16x16_t c6 = _mm256_unpackhi_epi32(b4,_zero);
    int16x16_t c5 = _mm256_unpacklo_epi32(b5,_zero);
    int16x16_t c7 = _mm256_unpackhi_epi32(b5,_zero);
    // int16x16_t a8 = _zero;
    // int16x16_t a9 = _zero;
    // int16x16_t b8 = _zero;
    // int16x16_t b9 = _zero;
    // int16x16_t a10 = _zero;
    // int16x16_t a11 = _zero;
    // int16x16_t b10 = _zero;
    // int16x16_t b11 = _zero;
    // int16x16_t c8 = _zero;
    // int16x16_t c10 = _zero;
    // int16x16_t c9 = _zero;
    // int16x16_t c11 = _zero;
    // int16x16_t a12 = _zero;
    // int16x16_t a13 = _zero;
    // int16x16_t b12 = _zero;
    // int16x16_t b13 = _zero;
    // int16x16_t a14 = _zero;
    // int16x16_t a15 = _zero;
    // int16x16_t b14 = _zero;
    // int16x16_t b15 = _zero;
    // int16x16_t c12 = _zero;
    // int16x16_t c14 = _zero;
    // int16x16_t c13 = _zero;
    // int16x16_t c15 = _zero;

    int16x16_t d0 = _mm256_unpacklo_epi64(c0,c4);
    int16x16_t d4 = _mm256_unpackhi_epi64(c0,c4);
    // int16x16_t d8 = _zero;
    // int16x16_t d12 = _zero;
    int16x16_t e0 = _mm256_permute2x128_si256(d0,_zero,0x20);
    int16x16_t e8 = _mm256_permute2x128_si256(d0,_zero,0x31);
    int16x16_t e4 = _mm256_permute2x128_si256(d4,_zero,0x20);
    int16x16_t e12 = _mm256_permute2x128_si256(d4,_zero,0x31);
    store_int16x16((int16x16_t*)(out[0]),  e0);
    store_int16x16((int16x16_t*)(out[8]),  e8);
    store_int16x16((int16x16_t*)(out[1]),  e4);
    store_int16x16((int16x16_t*)(out[9]), e12);

    int16x16_t d1 = _mm256_unpacklo_epi64(c1,c5);
    int16x16_t d5 = _mm256_unpackhi_epi64(c1,c5);
    // int16x16_t d9 = _zero;
    // int16x16_t d13 = _zero;
    int16x16_t e1 = _mm256_permute2x128_si256(d1,_zero,0x20);
    int16x16_t e9 = _mm256_permute2x128_si256(d1,_zero,0x31);
    int16x16_t e5 = _mm256_permute2x128_si256(d5,_zero,0x20);
    int16x16_t e13 = _mm256_permute2x128_si256(d5,_zero,0x31);
    store_int16x16((int16x16_t*)(out[ 4]),  e1);
    store_int16x16((int16x16_t*)(out[12]),  e9);
    store_int16x16((int16x16_t*)(out[ 5]),  e5);
    store_int16x16((int16x16_t*)(out[13]), e13);

    int16x16_t d2 = _mm256_unpacklo_epi64(c2,c6);
    int16x16_t d6 = _mm256_unpackhi_epi64(c2,c6);
    // int16x16_t d10 = _zero;
    // int16x16_t d14 = _zero;
    int16x16_t e2 = _mm256_permute2x128_si256(d2,_zero,0x20);
    int16x16_t e10 = _mm256_permute2x128_si256(d2,_zero,0x31);
    int16x16_t e6 = _mm256_permute2x128_si256(d6,_zero,0x20);
    int16x16_t e14 = _mm256_permute2x128_si256(d6,_zero,0x31);
    store_int16x16((int16x16_t*)(out[ 2]),  e2);
    store_int16x16((int16x16_t*)(out[ 3]),  e6);
    store_int16x16((int16x16_t*)(out[10]), e10);
    store_int16x16((int16x16_t*)(out[11]), e14);

    int16x16_t d3 = _mm256_unpacklo_epi64(c3,c7);
    int16x16_t d7 = _mm256_unpackhi_epi64(c3,c7);
    // int16x16_t d11 = _zero;
    // int16x16_t d15 = _zero;
    int16x16_t e3 = _mm256_permute2x128_si256(d3,_zero,0x20);
    int16x16_t e11 = _mm256_permute2x128_si256(d3,_zero,0x31);
    int16x16_t e7 = _mm256_permute2x128_si256(d7,_zero,0x20);
    int16x16_t e15 = _mm256_permute2x128_si256(d7,_zero,0x31);
    store_int16x16((int16x16_t*)(out[ 6]),  e3);
    store_int16x16((int16x16_t*)(out[ 7]),  e7);
    store_int16x16((int16x16_t*)(out[14]), e11);
    store_int16x16((int16x16_t*)(out[15]), e15);

}

void transpose_partial3_inv(int16_t *out0, int16_t *out1, int16_t *out2, int16_t *out3, int16_t *out4, int16_t *out5, const int16_t in[16][16])
{
    int16x16_t a0 = load_int16x16((int16x16_t*)&in[0][0]);
    int16x16_t a1 = load_int16x16((int16x16_t*)&in[1][0]);
    int16x16_t b0 = _mm256_unpacklo_epi16(a0,a1);
    int16x16_t b1 = _mm256_unpackhi_epi16(a0,a1);
    int16x16_t a2 = load_int16x16((int16x16_t*)&in[2][0]);
    int16x16_t a3 = load_int16x16((int16x16_t*)&in[3][0]);
    int16x16_t b2 = _mm256_unpacklo_epi16(a2,a3);
    int16x16_t b3 = _mm256_unpackhi_epi16(a2,a3);
    int16x16_t c0 = _mm256_unpacklo_epi32(b0,b2);
    int16x16_t c2 = _mm256_unpackhi_epi32(b0,b2);
    int16x16_t c1 = _mm256_unpacklo_epi32(b1,b3);
    // int16x16_t c3 = _mm256_unpackhi_epi32(b1,b3);
    int16x16_t a4 = load_int16x16((int16x16_t*)&in[4][0]);
    int16x16_t a5 = load_int16x16((int16x16_t*)&in[5][0]);
    int16x16_t b4 = _mm256_unpacklo_epi16(a4,a5);
    int16x16_t b5 = _mm256_unpackhi_epi16(a4,a5);
    int16x16_t a6 = load_int16x16((int16x16_t*)&in[6][0]);
    int16x16_t a7 = load_int16x16((int16x16_t*)&in[7][0]);
    int16x16_t b6 = _mm256_unpacklo_epi16(a6,a7);
    int16x16_t b7 = _mm256_unpackhi_epi16(a6,a7);
    int16x16_t c4 = _mm256_unpacklo_epi32(b4,b6);
    int16x16_t c6 = _mm256_unpackhi_epi32(b4,b6);
    int16x16_t c5 = _mm256_unpacklo_epi32(b5,b7);
    // int16x16_t c7 = _mm256_unpackhi_epi32(b5,b7);
    int16x16_t a8 = load_int16x16((int16x16_t*)&in[8][0]);
    int16x16_t a9 = load_int16x16((int16x16_t*)&in[9][0]);
    int16x16_t b8 = _mm256_unpacklo_epi16(a8,a9);
    int16x16_t b9 = _mm256_unpackhi_epi16(a8,a9);
    int16x16_t a10 = load_int16x16((int16x16_t*)&in[10][0]);
    int16x16_t a11 = load_int16x16((int16x16_t*)&in[11][0]);
    int16x16_t b10 = _mm256_unpacklo_epi16(a10,a11);
    int16x16_t b11 = _mm256_unpackhi_epi16(a10,a11);
    int16x16_t c8 = _mm256_unpacklo_epi32(b8,b10);
    int16x16_t c10 = _mm256_unpackhi_epi32(b8,b10);
    int16x16_t c9 = _mm256_unpacklo_epi32(b9,b11);
    // int16x16_t c11 = _mm256_unpackhi_epi32(b9,b11);
    int16x16_t a12 = load_int16x16((int16x16_t*)&in[12][0]);
    int16x16_t a13 = load_int16x16((int16x16_t*)&in[13][0]);
    int16x16_t b12 = _mm256_unpacklo_epi16(a12,a13);
    int16x16_t b13 = _mm256_unpackhi_epi16(a12,a13);
    int16x16_t a14 = load_int16x16((int16x16_t*)&in[14][0]);
    int16x16_t a15 = load_int16x16((int16x16_t*)&in[15][0]);
    int16x16_t b14 = _mm256_unpacklo_epi16(a14,a15);
    int16x16_t b15 = _mm256_unpackhi_epi16(a14,a15);
    int16x16_t c12 = _mm256_unpacklo_epi32(b12,b14);
    int16x16_t c14 = _mm256_unpackhi_epi32(b12,b14);
    int16x16_t c13 = _mm256_unpacklo_epi32(b13,b15);
    // int16x16_t c15 = _mm256_unpackhi_epi32(b13,b15);

    int16x16_t d0 = _mm256_unpacklo_epi64(c0,c4);
    int16x16_t d4 = _mm256_unpackhi_epi64(c0,c4);
    int16x16_t d8 = _mm256_unpacklo_epi64(c8,c12);
    int16x16_t d12 = _mm256_unpackhi_epi64(c8,c12);
    int16x16_t e0 = _mm256_permute2x128_si256(d0,d8,0x20);
    // int16x16_t e8 = _mm256_permute2x128_si256(d0,d8,0x31);
    int16x16_t e4 = _mm256_permute2x128_si256(d4,d12,0x20);
    // int16x16_t e12 = _mm256_permute2x128_si256(d4,d12,0x31);
    store_int16x16((int16x16_t*)(out0),  e0);
    // store_int16x16((int16x16_t*)(out[8]),  e8);
    store_int16x16((int16x16_t*)(out1),  e4);
    // store_int16x16((int16x16_t*)(out[9]), e12);

    int16x16_t d1 = _mm256_unpacklo_epi64(c1,c5);
    int16x16_t d5 = _mm256_unpackhi_epi64(c1,c5);
    int16x16_t d9 = _mm256_unpacklo_epi64(c9,c13);
    int16x16_t d13 = _mm256_unpackhi_epi64(c9,c13);
    int16x16_t e1 = _mm256_permute2x128_si256(d1,d9,0x20);
    // int16x16_t e9 = _mm256_permute2x128_si256(d1,d9,0x31);
    int16x16_t e5 = _mm256_permute2x128_si256(d5,d13,0x20);
    // int16x16_t e13 = _mm256_permute2x128_si256(d5,d13,0x31);
    store_int16x16((int16x16_t*)(out4),  e1);
    // store_int16x16((int16x16_t*)(out[12]),  e9);
    store_int16x16((int16x16_t*)(out5),  e5);
    // store_int16x16((int16x16_t*)(out[13]), e13);

    int16x16_t d2 = _mm256_unpacklo_epi64(c2,c6);
    int16x16_t d6 = _mm256_unpackhi_epi64(c2,c6);
    int16x16_t d10 = _mm256_unpacklo_epi64(c10,c14);
    int16x16_t d14 = _mm256_unpackhi_epi64(c10,c14);
    int16x16_t e2 = _mm256_permute2x128_si256(d2,d10,0x20);
    // int16x16_t e10 = _mm256_permute2x128_si256(d2,d10,0x31);
    int16x16_t e6 = _mm256_permute2x128_si256(d6,d14,0x20);
    // int16x16_t e14 = _mm256_permute2x128_si256(d6,d14,0x31);
    store_int16x16((int16x16_t*)(out2),  e2);
    store_int16x16((int16x16_t*)(out3),  e6);
    // store_int16x16((int16x16_t*)(out[10]), e10);
    // store_int16x16((int16x16_t*)(out[11]), e14);

    // int16x16_t d3 = _mm256_unpacklo_epi64(c3,c7);
    // int16x16_t d7 = _mm256_unpackhi_epi64(c3,c7);
    // int16x16_t d11 = _mm256_unpacklo_epi64(c11,c15);
    // int16x16_t d15 = _mm256_unpackhi_epi64(c11,c15);
    // int16x16_t e3 = _mm256_permute2x128_si256(d3,d11,0x20);
    // int16x16_t e11 = _mm256_permute2x128_si256(d3,d11,0x31);
    // int16x16_t e7 = _mm256_permute2x128_si256(d7,d15,0x20);
    // int16x16_t e15 = _mm256_permute2x128_si256(d7,d15,0x31);
    // store_int16x16((int16x16_t*)(out[ 6]),  e3);
    // store_int16x16((int16x16_t*)(out[ 7]),  e7);
    // store_int16x16((int16x16_t*)(out[14]), e11);
    // store_int16x16((int16x16_t*)(out[15]), e15);

}

