
# Compilation
```
make
```

# Correctness
Type
```
./test
```
Sample output:
```
polymul finished!
ntrup mul finished!
```

# Unit Tests
Type
```
./unit_test
```

Sample output:
```
================================
asm schoolbook 2x2 passed!
asm schoolbook cyclic 2x2 passed!
asm schoolbook negacyclic 2x2 passed!
asm FFT cyclic 2x2 passed!
================================
asm schoolbook 4x4 passed!
asm schoolbook cyclic 4x4 passed!
asm schoolbook negacyclic 4x4 passed!
asm FFT cyclic 4x4 passed!
asm FFT negacyclic 4x4 passed!
================================
asm Karatsuba 8x8 passed!
asm FFT cyclic 8x8 passed!
asm Karatsuba negacyclic 8x8 passed!
================================
asm FFT cyclic precompute 16x16 passed!
asm weighted Karatsuba and FFT2 are compatible (16x16)!
================================
```

# Benchmark
Type
```
./bench
```
Sample output:
```
polymul cycles: 13876, 13892, 13912
ntrup mul cycles: 14288, 14312, 14336
```

# Microbenchmark
Type
```
./microbench
```

Sample output:
```
8 repetitions
======== 2x2 ========
asm schoolbook 2x2: 120, 120, 120
asm schoolbook cyclic 2x2: 120, 120, 120
asm schoolbook negacyclic 2x2: 120, 120, 120
======== 4x4 ========
asm schoolbook 4x4: 476, 476, 480
asm schoolbook cyclic 4x4: 496, 500, 504
asm schoolbook negacyclic 4x4: 496, 500, 504
asm FFT cyclic 4x4: 256, 260, 260
asm FFT negacyclic 4x4: 612, 616, 620
======== 8x8 ========
asm Karatsuba 8x8: 1600, 1608, 1620
asm FFT cyclic 8x8: 944, 948, 948
asm Karatsuba negacyclic 8x8: 1644, 1648, 1656
======== 16x16 ========
asm FFT cyclic 16x16: 2508, 2512, 2516
asm weighted 16x16: 5564, 5584, 5604
asm weighted doubling 16x16: 5872, 5888, 5908
CT weighted 16x16: 4392, 4400, 4408
======== Butterflies ========
asm 3x2: 5992, 6028, 6052
asm 3x2 pre: 5020, 5024, 5032
asm 3x2 post: 3512, 3552, 3576
asm Rader-17: 2588, 2592, 2600
asm Rader-17 scaled indices: 2592, 2596, 2600
======== Transpose ========
transpose: 560, 560, 560
transpose partial 3: 364, 368, 368
transpose partial 3 inv: 412, 416, 416
```
