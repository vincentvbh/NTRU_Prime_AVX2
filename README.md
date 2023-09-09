
# Polynomial Multiplication of NTRU Prime with AVX2


This repository accompanies with the paper **Pushing the Limit of Vectorized Polynomial Multiplication for NTRU Prime**.
We target the big-by-big polynomial multiplication used in <tt>sntrup761</tt>.
Compared to the state-of-the-art AVX2-optimized implementation by [OpenSSLNTRU: Faster post-quantum TLS key exchange](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein),
Our big-by-big polynomial multiplication is 1.99 times faster on Haswell and 2.16 times faster on Skylake.


For the batch key generation with batch size 32, the amortized cost of key generation is reduced by $12\%$ on Haswell and $8 \%$ on Skylake.
For the encapsulation, the performance cycles are reduced by $7 \%$ on Haswell and $10 \%$ on Skylake.
As for the decapsulation, the performance cycles are reduced by $10 \%$ on Haswell and $13 \%$ on Skylake.

## Performance of Polynomial Multiplications

### Haswell Results

| Operation | [BBCT22](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein) (our benchmark) | This work |
| ------ | -------------- | ----------------------- |
| <tt>mulcore</tt> | 23,460 | 12,336 |
| <tt>polymul</tt> | 25,356 | 12,760 |

### Skylake Results

| Operation | [BBCT22](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein) (our benchmark) | This work |
| ------ | -------------- | ----------------------- |
| <tt>mulcore</tt> | 20,070 | 9,778 |
| <tt>polymul</tt> | 21,364 | 9,876 |


## Performance of Scheme <tt>sntrup761</tt>

### Haswell Results

| Operation | [BBCT22](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein) (our benchmark)  | This work |
| ------ | -------------- | ----------------------- |
| Batch key generation (batch size 32, package <tt>libsntrup761-20210608</tt>) | 154,552 | 136,003 |
| Encapsulation (package <tt>supercop-20230530</tt>) | 47,464 | 44,108 |
| Decapsulation (package <tt>supercop-20230530</tt>) | 56,064 | 50,080 |

### Skylake Results

| Operation | [BBCT22](https://www.usenix.org/conference/usenixsecurity22/presentation/bernstein) (our benchmark)  | This work |
| ------ | -------------- | ----------------------- |
| Batch key generation (batch size 32, package <tt>libsntrup761-20210608</tt>) | 129,159 | 118,939 |
| Encapsulation (package <tt>supercop-20230530</tt>) | 40,653 | 36,486 |
| Decapsulation (package <tt>supercop-20230530</tt>) | 47,387 | 41,070|

## Benchmarking Platform

### Haswell
```
Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         39 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  4
  On-line CPU(s) list:   0-3
Vendor ID:               GenuineIntel
  Model name:            Intel(R) Core(TM) i7-4770K CPU @ 3.50GHz
    CPU family:          6
    Model:               60
    Thread(s) per core:  1
    Core(s) per socket:  4
    Socket(s):           1
    Stepping:            3
    BogoMIPS:            6983.29
    Flags:               fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb r
                         dtscp lm constant_tsc arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc cpuid aperfmperf pni pclmulqdq dtes64 monitor ds_cpl vmx est
                          tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid sse4_1 sse4_2 movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm cpuid_fault epb
                         invpcid_single pti ssbd ibrs ibpb stibp tpr_shadow vnmi flexpriority ept vpid ept_ad fsgsbase tsc_adjust bmi1 avx2 smep bmi2 erms invpcid xs
                         aveopt dtherm arat pln pts md_clear flush_l1d
Virtualization features:
  Virtualization:        VT-x
Caches (sum of all):
  L1d:                   128 KiB (4 instances)
  L1i:                   128 KiB (4 instances)
  L2:                    1 MiB (4 instances)
  L3:                    8 MiB (1 instance)
NUMA:
  NUMA node(s):          1
  NUMA node0 CPU(s):     0-3
Vulnerabilities:
  Itlb multihit:         KVM: Mitigation: VMX disabled
  L1tf:                  Mitigation; PTE Inversion; VMX conditional cache flushes, SMT disabled
  Mds:                   Mitigation; Clear CPU buffers; SMT disabled
  Meltdown:              Mitigation; PTI
  Mmio stale data:       Unknown: No mitigations
  Retbleed:              Not affected
  Spec store bypass:     Mitigation; Speculative Store Bypass disabled via prctl and seccomp
  Spectre v1:            Mitigation; usercopy/swapgs barriers and __user pointer sanitization
  Spectre v2:            Mitigation; Retpolines, IBPB conditional, IBRS_FW, RSB filling, PBRSB-eIBRS Not affected
  Srbds:                 Mitigation; Microcode
  Tsx async abort:       Not affected
```

### Skylake
```
Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         39 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  8
  On-line CPU(s) list:   0-7
Vendor ID:               GenuineIntel
  Model name:            Intel(R) Xeon(R) CPU E3-1275 v5 @ 3.60GHz
    CPU family:          6
    Model:               94
    Thread(s) per core:  2
    Core(s) per socket:  4
    Socket(s):           1
    Stepping:            3
    CPU max MHz:         4000.0000
    CPU min MHz:         800.0000
    BogoMIPS:            7200.00
    Flags:               fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mc
                         a cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss
                         ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc art
                          arch_perfmon pebs bts rep_good nopl xtopology nonstop_
                         tsc cpuid aperfmperf pni pclmulqdq dtes64 monitor ds_cp
                         l vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid ss
                         e4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes
                         xsave avx f16c rdrand lahf_lm abm 3dnowprefetch cpuid_f
                         ault epb invpcid_single pti ssbd ibrs ibpb stibp tpr_sh
                         adow vnmi flexpriority ept vpid ept_ad fsgsbase tsc_adj
                         ust bmi1 avx2 smep bmi2 erms invpcid mpx rdseed adx sma
                         p clflushopt intel_pt xsaveopt xsavec xgetbv1 xsaves dt
                         herm ida arat pln pts hwp hwp_notify hwp_act_window hwp
                         _epp md_clear flush_l1d arch_capabilities
Virtualization features:
  Virtualization:        VT-x
Caches (sum of all):
  L1d:                   128 KiB (4 instances)
  L1i:                   128 KiB (4 instances)
  L2:                    1 MiB (4 instances)
  L3:                    8 MiB (1 instance)
NUMA:
  NUMA node(s):          1
  NUMA node0 CPU(s):     0-7
Vulnerabilities:
  Itlb multihit:         KVM: Mitigation: VMX disabled
  L1tf:                  Mitigation; PTE Inversion; VMX conditional cache flushe
                         s, SMT vulnerable
  Mds:                   Mitigation; Clear CPU buffers; SMT vulnerable
  Meltdown:              Mitigation; PTI
  Mmio stale data:       Mitigation; Clear CPU buffers; SMT vulnerable
  Retbleed:              Mitigation; IBRS
  Spec store bypass:     Mitigation; Speculative Store Bypass disabled via prctl
                          and seccomp
  Spectre v1:            Mitigation; usercopy/swapgs barriers and __user pointer
                          sanitization
  Spectre v2:            Mitigation; IBRS, IBPB conditional, STIBP conditional,
                         RSB filling, PBRSB-eIBRS Not affected
  Srbds:                 Mitigation; Microcode
  Tsx async abort:       Mitigation; TSX disabled
```
