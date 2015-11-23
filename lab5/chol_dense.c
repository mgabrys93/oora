#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


/* for Fortran - j*n + i */
//#define IDX(i, j, n)	(((i) * (n)) + (j))
#define IDX(i, j, n) (((j)+ (i)*(n)))
#define BLOCK_SIZE 8
#define SIZE 800

static double gtod_ref_time_sec = 0.0;

int
chol(double *A, unsigned int n)
{
	unsigned int i;
	unsigned int j;
	unsigned int k;

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			for (k = 0; k < j; ++k) {
				A[IDX(i, j, n)] -= A[IDX(i, k, n)] *
				    A[IDX(j, k, n)];
			}
		}

		if (A[IDX(j, j, n)] < 0.0) {
			return (1);
		}

		A[IDX(j, j, n)] = sqrt(A[IDX(j, j, n)]);
		for (i = j + 1; i < n; i++)
			A[IDX(i, j, n)] /= A[IDX(j, j, n)];
	}

	return (0);
}

int mm(double first[SIZE], double second[SIZE], double multiply[SIZE])
{
  int i,j,k; 
  double sum = 0;
  for (i = 0; i < SIZE; i++) { //rows in multiply
    for (j = 0; j < SIZE; j++) { //columns in multiply
      for (k = 0; k < SIZE; k++) { //columns in first and rows in second
	    sum = sum + first[IDX(i, k, SIZE)]*second[IDX(k, j, SIZE)];
	  } 
          multiply[IDX(i, j, SIZE)] = sum;
	  sum = 0;
    }
  }
  return 0;
}

double *generateSPDmatrix(){
	double *A = malloc(SIZE * SIZE * sizeof(double));
	double *At = malloc(SIZE * SIZE * sizeof(double));
	double *SPDM = malloc(SIZE * SIZE * sizeof(double));
	
	unsigned i, j;

	for(i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			double tmp = ((double)rand())/RAND_MAX;
			A[IDX(i, j, SIZE)] = tmp;
			A[IDX(j, i, SIZE)] = tmp;
		}
	}

	mm(A, At, SPDM);

	return SPDM;
}

double dclock()
{
  double         the_time, norm_sec;
  struct timeval tv;
  gettimeofday( &tv, NULL );
  if ( gtod_ref_time_sec == 0.0 )
    gtod_ref_time_sec = ( double ) tv.tv_sec;
  norm_sec = ( double ) tv.tv_sec - gtod_ref_time_sec;
  the_time = norm_sec + tv.tv_usec * 1.0e-6;
  return the_time;
}

int
main()
{
	srand((unsigned int)time(NULL));
	double *A;
	double dtime;
	int i, j, n, ret;
	
	n = 3;
	A = generateSPDmatrix();

	dtime = dclock();
 	chol(A, SIZE);
 	dtime = dclock()-dtime;
 	double gflops = ((1.0/3.0) * SIZE * SIZE * SIZE * 10e-9) / dtime;
	printf( "Time: %le \n", dtime);
	printf("Gflops: %le \n", gflops);
	return 0;
}
