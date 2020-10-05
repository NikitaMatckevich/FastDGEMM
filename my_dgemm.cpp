#include <cstdlib>
#include <cstring>
#include <string>
#include <charconv>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <immintrin.h>
using namespace std;
using thread_storage_t = 
  typename aligned_storage<sizeof(thread),alignof(thread)>::type;

constexpr int UNROLL = 4;
constexpr int BLOCKSIZE = 32;
constexpr int MAX_NUM_THREADS = 4;

void do_block(const int N,
              const double* A, const double* B, double* C)
{
  for (int j = 0; j < BLOCKSIZE; j += 4*UNROLL) {
    for (int i = 0; i < BLOCKSIZE; i++) {
      __m256d c[UNROLL];
      for (int x = 0; x < UNROLL; x++)
        c[x] = _mm256_load_pd(C + i*N + j + x*4);
      for (int k = 0; k < BLOCKSIZE; k++) {
        __m256d a = _mm256_broadcast_sd(A + i*N + k);
        for (int x = 0; x < UNROLL; x++)
        c[x] = _mm256_add_pd(c[x],
            _mm256_mul_pd(_mm256_load_pd(B + k*N + j + x*4), a));
      }
      for (int x = 0; x < UNROLL; x++)
        _mm256_store_pd(C + i*N + j + x*4, c[x]);
    }
  }
}

void loc_my_dgemm(const int M, const int N,
                  const double* A, const double* B, double* C)
{
  for (int si = 0; si < M; si += BLOCKSIZE)
    for (int sj = 0; sj < N; sj += BLOCKSIZE)
      for (int sk = 0; sk < N; sk += BLOCKSIZE)
        do_block(N, A + N*si + sk, B + N*sk + sj, C + N*si + sj);
}

void glob_my_dgemm(const int NUM_THREADS, const int N,
                   const double* A, const double* B, double* C)
{
  thread_storage_t thread_pool[MAX_NUM_THREADS];
  int wN = N / NUM_THREADS;
  for (int th_id = 0; th_id < NUM_THREADS-1; th_id++) {
    int offset = (th_id + 1) * wN * N;
    auto wA = A + offset;
    auto wC = C + offset;
    new (&thread_pool[th_id]) thread(loc_my_dgemm, wN, N, wA, B, wC);
  }
  loc_my_dgemm(wN, N, A, B, C);
  for (int th_id = 0; th_id < NUM_THREADS-1; th_id++)
    reinterpret_cast<thread*>(&thread_pool[th_id])->join();
}

void run(const int NUM_THREADS, const int N,
         const double* A, const double* B, double* C,
         double* duration, double* gflops)
{
  constexpr int nbTests = 10;
  const double iNbTests = 1./nbTests;

  *duration = 0.;
  *gflops   = 0.;

  const int N2 = N*N;

  for (int test = 0; test < nbTests; test++) {

    auto t_beg = chrono::high_resolution_clock::now();
    glob_my_dgemm(NUM_THREADS, N, A, B, C);
    auto t_end = chrono::high_resolution_clock::now();

    double test_d = chrono::duration<double,std::milli>(t_end - t_beg).count();
    double test_g = 2. * N * N * N;
    test_g = test_g/test_d*1.e-6;
    
    cout << test << ":\t" << test_d << "ms\t" << test_g << "gflop/s\n";
    /*  // OUTPUT:
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++)
        cout << C[i*N + j] << " ";
      cout << "\n";
    }*/
    
    memset(C, 0, sizeof(double)*N2);

    *duration += iNbTests*test_d;
    *gflops   += iNbTests*test_g;
  }
}

int main(int argc, char** argv) {
try {
  if (argc != 3)
    throw invalid_argument("invalid number of arguments: expected 2");
  int NUM_THREADS;
  {
    string str(argv[1]);
    auto [p, ec] = from_chars(str.data(), str.data() + str.size(), NUM_THREADS);
    if (ec != std::errc())
      throw invalid_argument("couldn't parse number of threads");
  }
  int N;
  {
    string str(argv[2]);
    auto [p, ec] = from_chars(str.data(), str.data() + str.size(), N);
    if (ec != std::errc())
      throw invalid_argument("couldn't parse number of matrix rows");
  }
  int N2 = N*N;

  // memory allocation for the arrays (alignment for m256 instructions)
  double *A = static_cast<double*>(aligned_alloc(alignof(double)*4, sizeof(double)*N2));
  double *B = static_cast<double*>(aligned_alloc(alignof(double)*4, sizeof(double)*N2));
  double *C = static_cast<double*>(aligned_alloc(alignof(double)*4, sizeof(double)*N2));
  // set arrays equal to zero
  memset(A, 0, sizeof(double)*N2);
  memset(B, 0, sizeof(double)*N2);
  memset(C, 0, sizeof(double)*N2);
  // diagonal test
  for (int i = 0; i < N2; i++) {
    double elem = (double)(i+1);
    A[i] =  elem;
    B[i] = -elem;
  }
    
  double duration = 0.;
  double gflops = 0.;
  run(NUM_THREADS, N, A, B, C, &duration, &gflops);
  
  ofstream fout;
  fout.open("myDGEMM.txt", std::ios_base::app);
  fout << N << "\t" << duration << "ms\t" << gflops << "gflop/s\n";
  fout.close();

  std::free(A);
  std::free(B);
  std::free(C);
  return 0;
}
catch (const exception& e) {
  cerr << e.what() << "\n";
}
}
