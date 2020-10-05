# FastDGEMM
Simple multithread version of Hennessy-Patterson optimized dgemm routine on C++.

Short benchmark code for the OpenBLAS library is added for comparison.

I would be happy to know your opinion about this experiment and get some feedback about the general methodology, benchmark reliability, and quality of the code itself.

I would also appreciate any idea on how to improve the performance of this code.

Stages of improvement ("Getting faster");

![alt](https://github.com/NikitaMatckevich/FastDGEMM/benchmarks/1-Compiler.png)

> Compiler optimizations

![alt](https://github.com/NikitaMatckevich/FastDGEMM/benchmarks/2-SIMD.png)

> AVX x86 intrinsics

![alt](https://github.com/NikitaMatckevich/FastDGEMM/benchmarks/3-Pipelined.png)

> Loop unrolling (-O3 compiler option shoul be enabled)

![alt](https://github.com/NikitaMatckevich/FastDGEMM/benchmarks/4-Multithread.png)

> New feature, realization maybe isn't perfect

![alt](https://github.com/NikitaMatckevich/FastDGEMM/benchmarks/5-CacheBlocking.png)

> Cache blocking (controlling the age of the array accesses)

![alt](https://github.com/NikitaMatckevich/FastDGEMM/benchmarks/6-OpenBLAS.png)

> Comparison with OpenBLAS. We've got a really big trip ahead...