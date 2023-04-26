#ifndef __AVX2_BASEMUL_H
#define __AVX2_BASEMUL_H

#include <stdint.h>

#include "params.h"
#include "__avx2_wrap.h"
#include "NTT_params.h"

extern int16_t const_buff[80];

extern void __asm_schoolbook2(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);
extern void __asm_cyclic_schoolbook2(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);
extern void __asm_negacyclic_schoolbook2(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);
extern void __asm_schoolbook4(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);

extern void __asm_cyclic_schoolbook4(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);
extern void __asm_negacyclic_schoolbook4(int16_t *des, int16_t *src1, int16_t *src2, int16_t *_const_buff);

void __barrett_int16x16(int16_t *a);

static inline
int16x16_t montmulmod_int16x16(int16x16_t src1, int16x16_t src2){


    int16x16_t lo, hi;
    int16x16_t qx16, qinvx16;

    qx16 = dup_const_int16x16(Q);
    qinvx16 = dup_const_int16x16(Qinv);

    lo = mullo_int16x16(src2, qinvx16);
    hi = mulhi_int16x16(src1, src2);
    lo = mullo_int16x16(src1, lo);
    lo = mulhi_int16x16(lo, qx16);

    return sub_int16x16(hi, lo);
}

static inline
int16x16_t montmulmod_precompute_int16x16(int16x16_t src1, int16x16_t src2, int16x16_t src2inv){


    int16x16_t lo, hi;
    int16x16_t qx16;

    qx16 = dup_const_int16x16(Q);

    hi = mulhi_int16x16(src1, src2);
    lo = mullo_int16x16(src1, src2inv);
    lo = mulhi_int16x16(lo, qx16);

    return sub_int16x16(hi, lo);
}

static inline
int16x16_t barrett_int16x16(int16x16_t a)
{
  return sub_int16x16(a,mullo_int16x16(mulhrs_int16x16(a,dup_const_int16x16(Qbar)),dup_const_int16x16(Q)));
}

static inline
void schoolbook2_Rmod_int16x16(int16x16_t des[3], int16x16_t src1[2], int16x16_t src2[2]){

    des[1] = add_int16x16(montmulmod_int16x16(src1[0], src2[1]), montmulmod_int16x16(src1[1], src2[0]));
    des[0] = montmulmod_int16x16(src1[0], src2[0]);
    des[2] = montmulmod_int16x16(src1[1], src2[1]);

}

static inline
void cyclic_schoolbook2_Rmod_int16x16(int16x16_t des[2], int16x16_t src1[2], int16x16_t src2[2]){

    des[0] = add_int16x16(montmulmod_int16x16(src1[0], src2[0]), montmulmod_int16x16(src1[1], src2[1]));
    des[1] = add_int16x16(montmulmod_int16x16(src1[0], src2[1]), montmulmod_int16x16(src1[1], src2[0]));

}

static inline
void negacyclic_schoolbook2_Rmod_int16x16(int16x16_t des[2], int16x16_t src1[2], int16x16_t src2[2]){

    des[0] = sub_int16x16(montmulmod_int16x16(src1[0], src2[0]), montmulmod_int16x16(src1[1], src2[1]));
    des[1] = add_int16x16(montmulmod_int16x16(src1[0], src2[1]), montmulmod_int16x16(src1[1], src2[0]));

}

static inline
void negacyclic_schoolbook2_precompute_Rmod_int16x16(
    int16x16_t des[2], int16x16_t src1[2], int16x16_t src2[2], int16x16_t src2inv[2]){

    des[0] = sub_int16x16(montmulmod_precompute_int16x16(src1[0], src2[0], src2inv[0]),
                          montmulmod_precompute_int16x16(src1[1], src2[1], src2inv[1]));
    des[1] = add_int16x16(montmulmod_precompute_int16x16(src1[0], src2[1], src2inv[1]),
                          montmulmod_precompute_int16x16(src1[1], src2[0], src2inv[0]));

}

static inline
void schoolbook4_Rmod_int16x16(int16x16_t des[7], int16x16_t src1[4], int16x16_t src2[4]){

    des[3] = add_int16x16( add_int16x16(montmulmod_int16x16(src1[0], src2[3]), montmulmod_int16x16(src1[1], src2[2])),
                           add_int16x16(montmulmod_int16x16(src1[2], src2[1]), montmulmod_int16x16(src1[3], src2[0])) );
    des[2] = add_int16x16( montmulmod_int16x16(src1[0], src2[2]),
                           add_int16x16(montmulmod_int16x16(src1[1], src2[1]), montmulmod_int16x16(src1[2], src2[0])) );
    des[4] = add_int16x16( montmulmod_int16x16(src1[1], src2[3]),
                           add_int16x16(montmulmod_int16x16(src1[2], src2[2]), montmulmod_int16x16(src1[3], src2[1])) );
    des[1] = add_int16x16(montmulmod_int16x16(src1[0], src2[1]), montmulmod_int16x16(src1[1], src2[0]));
    des[5] = add_int16x16(montmulmod_int16x16(src1[2], src2[3]), montmulmod_int16x16(src1[3], src2[2]));
    des[0] = montmulmod_int16x16(src1[0], src2[0]);
    des[6] = montmulmod_int16x16(src1[3], src2[3]);
    des[3] = barrett_int16x16(des[3]);

}

static inline
void schoolbook4_precompute_Rmod_int16x16(
    int16x16_t des[7], int16x16_t src1[4], int16x16_t src2[4], int16x16_t src2inv[4]){

    des[3] = add_int16x16( add_int16x16(
                                montmulmod_int16x16(src1[0], src2[3]),
                                montmulmod_int16x16(src1[1], src2[2])),
                           add_int16x16(
                                montmulmod_int16x16(src1[2], src2[1]),
                                montmulmod_int16x16(src1[3], src2[0])) );
    des[2] = add_int16x16(
                                montmulmod_int16x16(src1[0], src2[2]),
                           add_int16x16(
                                montmulmod_int16x16(src1[1], src2[1]),
                                montmulmod_int16x16(src1[2], src2[0])) );
    des[4] = add_int16x16(
                                montmulmod_int16x16(src1[1], src2[3]),
                           add_int16x16(
                                montmulmod_int16x16(src1[2], src2[2]),
                                montmulmod_int16x16(src1[3], src2[1])) );
    des[1] = add_int16x16(
                                montmulmod_int16x16(src1[0], src2[1]),
                                montmulmod_int16x16(src1[1], src2[0]));
    des[5] = add_int16x16(
                                montmulmod_int16x16(src1[2], src2[3]),
                                montmulmod_int16x16(src1[3], src2[2]));
    des[0] = montmulmod_precompute_int16x16(src1[0], src2[0], src2inv[0]);
    des[6] = montmulmod_precompute_int16x16(src1[3], src2[3], src2inv[3]);
    des[3] = barrett_int16x16(des[3]);

}

static inline
void cyclic_schoolbook4_Rmod_int16x16(int16x16_t des[4], int16x16_t src1[4], int16x16_t src2[4]){

    des[0] = barrett_int16x16(add_int16x16(
                                montmulmod_int16x16(src1[0], src2[0]),
                                add_int16x16(add_int16x16(
                                             montmulmod_int16x16(src1[1], src2[3]), montmulmod_int16x16(src1[2], src2[2])),
                                             montmulmod_int16x16(src1[3], src2[1])) ));

    des[1] = barrett_int16x16(add_int16x16(
                                add_int16x16(montmulmod_int16x16(src1[0], src2[1]), montmulmod_int16x16(src1[1], src2[0])),
                                add_int16x16(montmulmod_int16x16(src1[2], src2[3]), montmulmod_int16x16(src1[3], src2[2]))
                                             ));

    des[2] = barrett_int16x16(add_int16x16(
                                add_int16x16(add_int16x16(
                                    montmulmod_int16x16(src1[0], src2[2]), montmulmod_int16x16(src1[1], src2[1])),
                                    montmulmod_int16x16(src1[2], src2[0])),
                                montmulmod_int16x16(src1[3], src2[3])));
    des[3] = barrett_int16x16(add_int16x16(
                                add_int16x16(montmulmod_int16x16(src1[0], src2[3]), montmulmod_int16x16(src1[1], src2[2])),
                                add_int16x16(montmulmod_int16x16(src1[2], src2[1]), montmulmod_int16x16(src1[3], src2[0]))
                                             ));

}

static inline
void negacyclic_schoolbook4_Rmod_int16x16(int16x16_t des[4], int16x16_t src1[4], int16x16_t src2[4]){

    des[0] = barrett_int16x16(sub_int16x16(
                                montmulmod_int16x16(src1[0], src2[0]),
                                add_int16x16(add_int16x16(
                                             montmulmod_int16x16(src1[1], src2[3]), montmulmod_int16x16(src1[2], src2[2])),
                                             montmulmod_int16x16(src1[3], src2[1])) ));

    des[1] = barrett_int16x16(sub_int16x16(
                                add_int16x16(montmulmod_int16x16(src1[0], src2[1]), montmulmod_int16x16(src1[1], src2[0])),
                                add_int16x16(montmulmod_int16x16(src1[2], src2[3]), montmulmod_int16x16(src1[3], src2[2]))
                                             ));

    des[2] = barrett_int16x16(sub_int16x16(
                                add_int16x16(add_int16x16(
                                    montmulmod_int16x16(src1[0], src2[2]), montmulmod_int16x16(src1[1], src2[1])),
                                    montmulmod_int16x16(src1[2], src2[0])),
                                montmulmod_int16x16(src1[3], src2[3])));
    des[3] = barrett_int16x16(add_int16x16(
                                add_int16x16(montmulmod_int16x16(src1[0], src2[3]), montmulmod_int16x16(src1[1], src2[2])),
                                add_int16x16(montmulmod_int16x16(src1[2], src2[1]), montmulmod_int16x16(src1[3], src2[0]))
                                             ));

}

static inline
void negacyclic_schoolbook4_precompute_Rmod_int16x16(
    int16x16_t des[4], int16x16_t src1[4], int16x16_t src2[4], int16x16_t src2inv[4]){

    des[0] = barrett_int16x16(sub_int16x16(
                                montmulmod_precompute_int16x16(src1[0], src2[0], src2inv[0]),
                                add_int16x16(add_int16x16(
                                             montmulmod_precompute_int16x16(src1[1], src2[3], src2inv[3]),
                                             montmulmod_precompute_int16x16(src1[2], src2[2], src2inv[2])),
                                             montmulmod_precompute_int16x16(src1[3], src2[1], src2inv[1])) ));

    des[1] = barrett_int16x16(sub_int16x16(
                                add_int16x16(montmulmod_precompute_int16x16(src1[0], src2[1], src2inv[1]),
                                             montmulmod_precompute_int16x16(src1[1], src2[0], src2inv[0])),
                                add_int16x16(montmulmod_precompute_int16x16(src1[2], src2[3], src2inv[3]),
                                             montmulmod_precompute_int16x16(src1[3], src2[2], src2inv[2]))
                                             ));

    des[2] = barrett_int16x16(sub_int16x16(
                                add_int16x16(add_int16x16(
                                    montmulmod_precompute_int16x16(src1[0], src2[2], src2inv[2]),
                                    montmulmod_precompute_int16x16(src1[1], src2[1], src2inv[1])),
                                    montmulmod_precompute_int16x16(src1[2], src2[0], src2inv[0])),
                                montmulmod_precompute_int16x16(src1[3], src2[3], src2inv[3])));
    des[3] = barrett_int16x16(add_int16x16(
                                add_int16x16(montmulmod_precompute_int16x16(src1[0], src2[3], src2inv[3]),
                                             montmulmod_precompute_int16x16(src1[1], src2[2], src2inv[2])),
                                add_int16x16(montmulmod_precompute_int16x16(src1[2], src2[1], src2inv[1]),
                                             montmulmod_precompute_int16x16(src1[3], src2[0], src2inv[0]))
                                             ));

}

#endif

