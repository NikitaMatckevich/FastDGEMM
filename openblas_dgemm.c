#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"
#include <cblas.h>

void run(const int N, const double* A, const double* B,
         double* C, double* duration, double* gflops)
{
  const int nbTests = 10;
  const double iNbTests = 1./nbTests;
  struct timeval start, finish;
  
  for (int test = 0; test < nbTests; test++) {
    srand((unsigned)time(NULL));
    printf ("%d run via CBLAS interface ...", test);

    gettimeofday(&start, NULL);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    N, N, N, 1.0, A, N, B, N, 1.0, C, N);
    gettimeofday(&finish, NULL);
    
    printf ("done.\n");

    double test_duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 +
                (double)(finish.tv_usec-start.tv_usec)) / 1000000;
    double test_gflops = 2. * N * N * N;
    test_gflops = test_gflops/test_duration*1.0e-9;
    printf("%d: \t %lfs \t %lfgflops\n\n", test, test_duration, test_gflops);

    *duration += iNbTests*test_duration;
    *gflops += iNbTests*test_gflops;
  }
}

int main(int argc, char** argv) {
  if (argc != 2) {
    printf("error");
    return 1;
  }
  int N = atoi(argv[1]);
  int N2 = N*N;
  //Memory allocation for the arrays:
  double *A,*B,*C;
  A = (double *)malloc( N2*sizeof( double ) );
  B = (double *)malloc( N2*sizeof( double ) );
  C = (double *)malloc( N2*sizeof( double ) );
  for (int i = 0; i < N2; i++) {
    A[i] = (double)(i+1);
    B[i] = (double)(-i-1);
    C[i] = 0.0;
  }
  
  double duration, gflops;
  run(N, A, B, C, &duration, &gflops);
  
  FILE *fp;
  fp = fopen("timeDGEMM.txt", "a");
  fprintf(fp, "%d \t %lfs \t %lfgflops\n", N, duration, gflops);
  fclose(fp);

  free(A);
  free(B);
  free(C);
  return 0;
}
