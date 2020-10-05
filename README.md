# FastDGEMM
Simple multithread version of Hennessy-Patterson optimized dgemm routine on C++.

Short benchmark code for the OpenBLAS library is added for comparison.

I would be happy to know your opinion about this experiment and get some feedback about the general methodology, benchmark reliability, and quality of the code itself.

I would also appreciate any idea on how to improve the performance of this code.

Stages of improvement ("Getting faster");

![](https://github.com/NikitaMatckevich/FastDGEMM/blob/master/benchmarks/1-Compiler.png)

> Compiler optimizations

![](https://github.com/NikitaMatckevich/FastDGEMM/blob/master/benchmarks/2-SIMD.png)

> AVX x86 intrinsics

![](https://github.com/NikitaMatckevich/FastDGEMM/blob/master/benchmarks/3-Pipelined.png)

> Loop unrolling (-O3 compiler option shoul be enabled)

![](https://github.com/NikitaMatckevich/FastDGEMM/blob/master/benchmarks/4-Multithread.png)

> New feature, realization maybe isn't perfect

![](https://github.com/NikitaMatckevich/FastDGEMM/blob/master/benchmarks/5-CacheBlocking.png)

> Cache blocking (controlling the age of the array accesses)

![](https://github.com/NikitaMatckevich/FastDGEMM/blob/master/benchmarks/6-OpenBLAS.png)

> Comparison with OpenBLAS. We've got a really big trip ahead...