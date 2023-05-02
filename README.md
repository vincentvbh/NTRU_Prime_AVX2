
# Polynomial Multiplication of NTRU Prime with AVX2

This repository accompanies with the paper [**Technical Report: Even Faster Polynomial Multiplication for NTRU Prime with AVX2**](https://eprint.iacr.org/2023/604).
Compared to the state-of-the-art optimized implementation by
[OpenSSLNTRU: Faster post-quantum TLS key exchange](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein),
our big-by-big polynomial multiplication for <tt>ntrulpr761/sntrup761</tt> is 1.77 times, 1.9 times, and 1.92 times faster on Haswell, Skylake, and Comet Lake.
Folder ``libsntrup`` contains their big-by-big polynomial multiplication and folder ``ntrup`` contains ours.
Please refer to the corresponding folders for instructions.

# Haswell Results

| Operation | [BBCT22](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein) | This work |
| ------ | -------------- | ----------------------- |
| <tt>mulcore</tt> | 23,460 | 13,892 |
| <tt>polymul</tt> | 25,356 | 14,312 |

# Skylake Results

| Operation | [BBCT22](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein) | This work |
| ------ | -------------- | ----------------------- |
| <tt>mulcore</tt> | 21,402 | 11,682 |
| <tt>polymul</tt> | 23,306 | 12,242 |

# Comet Lake Results

| Operation | [BBCT22](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein) | This work |
| ------ | -------------- | ----------------------- |
| <tt>mulcore</tt> | 16,154 | 8,570 |
| <tt>polymul</tt> | 16,852 | 8,776 |

