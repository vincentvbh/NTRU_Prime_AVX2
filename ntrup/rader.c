
#include "rader.h"

#include "naive_mult.h"
#include "ring.h"

#include "NTT_params.h"

#include "__avx2.h"

// g = 3
// map x to g^x
size_t rader_in_permute[17] = {
0,
3, 9, 10, 13,
5, 15, 11, 16,
14, 8, 7, 4,
12, 2, 6, 1
};
// map x to g^{16 - x}
size_t rader_out_permute[17] = {
0,
6, 2, 12, 4,
7, 8, 14, 16,
11, 15, 5, 13,
10, 9, 3, 1
};

void rader_17(int16_t *des, int16_t *src, int16_t *twiddle_table, size_t jump){

    int16_t src_buff[17];
    int16_t twiddle_table_buff[17];
    int16_t twiddle;
    int16_t head;

    for(size_t i = 1; i < 17; i++){
        src_buff[i] = src[rader_in_permute[i] * jump];
        twiddle_table_buff[i] = twiddle_table[rader_in_permute[17 - i]];
    }

    head = src[0];
    src_buff[0] = head;

    for(size_t i = 1; i < 17; i++){
        coeff_ring.addZ(src_buff, src_buff, src_buff + i);
    }

    twiddle = 1;
    naive_mulR(src_buff + 1, src_buff + 1, twiddle_table_buff + 1, 16, &twiddle, coeff_ring);

    for(size_t i = 1; i < 17; i++){
        coeff_ring.addZ(src_buff + i, src_buff + i, &head);
    }

    for(size_t i = 0; i < 17; i++){
        des[rader_out_permute[i] * jump] = src_buff[i];
    }


}













