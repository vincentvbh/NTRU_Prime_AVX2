
CC = gcc

CFLAGS += -mcpu=native -mtune=native -march=native -O3 -mavx2 -Wall -z noexecstack

COMMON_PATH = ../common

CFLAGS += -I$(COMMON_PATH)

COMMON_SOURCE = $(COMMON_PATH)/tools.c $(COMMON_PATH)/naive_mult.c $(COMMON_PATH)/gen_table.c $(COMMON_PATH)/ntt_c.c \
	$(COMMON_PATH)/cpucycles.c

ASM_HEADERs =
ASM_SOURCEs =
C_HEADERs = params.h NTT_params.h ring.h montproduct.h avx.h
C_SOURCEs = ring.c avx.c
AVX2_HEADERs =
AVX2_SOURCEs =

HEADERs = $(ASM_HEADERs) $(C_HEADERs) $(AVX2_HEADERs)
SOURCEs = $(ASM_SOURCEs) $(COMMON_SOURCE) $(C_SOURCEs) $(AVX2_SOURCEs)

all: test bench

test: $(addsuffix .c, $@) $(HEADERs) $(SOURCEs)
	$(CC) $(CFLAGS) $(SOURCEs) -o $@ $(addsuffix .c, $@)

bench: $(addsuffix .c, $@) $(HEADERs) $(SOURCEs)
	$(CC) $(CFLAGS) $(SOURCEs) -o $@ $(addsuffix .c, $@)


.PHONY: clean
clean:
	rm -f test
	rm -f bench
