#ifndef RADER_H
#define RADER_H

#include <stddef.h>
#include <stdint.h>

// g = 3
// map x to g^x
extern size_t rader_in_permute[17];
// map x to g^{16 - x}
extern size_t rader_out_permute[17];

void rader_17(int16_t *des, int16_t *src, int16_t *twiddle_table, size_t jump);



#endif

