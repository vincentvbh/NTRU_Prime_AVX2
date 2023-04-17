#ifndef __AVX2_FFT_H
#define __AVX2_FFT_H

#include "__avx2_wrap.h"
#include "__avx2_basemul.h"

extern int16_t const_twiddle51p[48];
extern int16_t const_twiddle51n[48];
extern int16_t const_twiddle102_last[8];

extern int16_t OMEGA3_buff[16];
extern int16_t OMEGA3INV_buff[16];

extern int16_t const_twiddle17[17 * 16];
extern int16_t const_twiddle17pre[17 * 16];

extern int16_t const_twiddle17inv[17 * 16];
extern int16_t const_twiddle17invpre[17 * 16];

extern uint64_t _tbl[17];
extern uint64_t permute_tbl[2 * 17 * 3];

extern int16_t twiddle_CT[48];
extern int16_t twiddle_GS[48];

extern void __asm_3x2(int16_t *srclo, int16_t *srchi, int16_t *_twiddle, int16_t *const_buff);
extern void __asm_3x2_pre(int16_t *srclo, int16_t *srchi, int16_t *_twiddle, int16_t *const_buff);
extern void __asm_3x2_post(int16_t *srclo, int16_t *srchi, int16_t *_twiddle, int16_t *const_buff);

extern void __asm_rader17(int16_t *des, const int16_t *src, int16_t *twiddle_table, int16_t *twiddle_inv_table, int16_t *_const_buff);
extern void __asm_rader17_scalei(int16_t *des, const int16_t *src, int16_t *twiddle_table, int16_t *twiddle_inv_table, int16_t *_const_buff);

extern void __asm_rader17_pre_tbl(int16_t *des, const int16_t *src, int16_t *twiddle_table, int16_t *twiddle_inv_table,
                                  int16_t *_const_buff, uint64_t *_tbl);
extern void __asm_rader17_post_tbl(int16_t *des, const int16_t *src, int16_t *twiddle_table, int16_t *twiddle_inv_table,
                                   int16_t *_const_buff, uint64_t *_tbl);

void rader_17_Rmod_int16x16(int16_t *des, int16_t *src, int16_t *twiddle_table, int16_t *twiddle_inv_table, size_t jump);
void rader_17_scalei_Rmod_int16x16(int16_t *des, int16_t *src, int16_t *twiddle_table, int16_t *twiddle_inv_table);

#endif

