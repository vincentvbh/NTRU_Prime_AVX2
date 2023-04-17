
CC = gcc

CFLAGS += -mcpu=native -mtune=native -march=native -O3 -mavx2 -Wall -z noexecstack

COMMON_PATH = ../common

CFLAGS += -I$(COMMON_PATH)

COMMON_SOURCE = $(COMMON_PATH)/tools.c $(COMMON_PATH)/naive_mult.c $(COMMON_PATH)/gen_table.c $(COMMON_PATH)/ntt_c.c \
	$(COMMON_PATH)/cpucycles.c

ASM_HEADERs = basemul_core.inc permute.inc
ASM_SOURCEs = basemul.S NTT.S
C_HEADERs = params.h NTT_params.h ring.h vec.h  rader.h
C_SOURCEs = ring.c vec.c  rader.c
AVX2_HEADERs = __avx2_wrap.h \
			   __avx2_permute.h __avx2_basemul.h __avx2_basemul_Karatsuba.h __avx2_basemul_FFT.h __avx2_FFT.h __avx2.h
AVX2_SOURCEs = __avx2_permute.c __avx2_basemul.c __avx2_basemul_Karatsuba.c __avx2_basemul_FFT.c __avx2_FFT.c __avx2.c

HEADERs = $(ASM_HEADERs) $(C_HEADERs) $(AVX2_HEADERs)
SOURCEs = $(ASM_SOURCEs) $(COMMON_SOURCE) $(C_SOURCEs) $(AVX2_SOURCEs)

all: test bench unit_test microbench

test: $(addsuffix .c, $@) $(HEADERs) $(SOURCEs)
	$(CC) $(CFLAGS) $(SOURCEs) -o $@ $(addsuffix .c, $@)

gen: $(addsuffix .c, $@) $(HEADERs) $(SOURCEs)
	$(CC) $(CFLAGS) $(SOURCEs) -o $@ $(addsuffix .c, $@)

unit_test: $(addsuffix .c, $@) $(HEADERs) $(SOURCEs)
	$(CC) $(CFLAGS) $(SOURCEs) -o $@ $(addsuffix .c, $@)

bench: $(addsuffix .c, $@) $(HEADERs) $(SOURCEs)
	$(CC) $(CFLAGS) $(SOURCEs) -o $@ $(addsuffix .c, $@)

microbench: $(addsuffix .c, $@) $(HEADERs) $(SOURCEs)
	$(CC) $(CFLAGS) $(SOURCEs) -o $@ $(addsuffix .c, $@)

.PHONY: clean
clean:
	rm -f test
	rm -f unit_test
	rm -f bench
	rm -f microbench
