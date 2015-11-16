#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "papi_fun.h"

/* for Fortran - j*n + i */
//#define IDX(i, j, n)	(((i) * (n)) + (j))
#define IDX(i, j, n) (((j)+ (i)*(n)))
#define BLOCK_SIZE 8
#define SIZE 800

static double gtod_ref_time_sec = 0.0;

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



int
chol_right_looking(double *A, unsigned int n)
{
	register unsigned int i, j, k;
	register double tmp = 0.0;

	for(k = 0; k < n; k++){
		A[IDX(k, k, n)] = sqrt(A[IDX(k, k, n)]);

		for(i = k+1; i < n; i++){
			tmp = A[IDX(k, k, n)];
			if(i + BLOCK_SIZE < n){
				A[IDX(i, k, n)] /= tmp;
				A[IDX(i+1, k, n)] /= tmp;
				A[IDX(i+2, k, n)] /= tmp;
				A[IDX(i+3, k, n)] /= tmp;
				A[IDX(i+4, k, n)] /= tmp;
				A[IDX(i+5, k, n)] /= tmp;
				A[IDX(i+6, k, n)] /= tmp;
				A[IDX(i+7, k, n)] /= tmp;
				i += 8;
			}else{
				A[IDX(i, k, n)] /= tmp;
				i++;
			}
		}

		for(j = k+1; j < n; j++){
			for(i = j; i < n; i++){
				A[IDX(i, j, n)] -= A[IDX(i, k, n)] * A[IDX(j, k, n)];
			}
		}
	}

	return (0);
}


int
chol_left_looking(double *A, unsigned int n)
{
	register unsigned int i, j, k;
	register double sum, tmp;
	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0.0;
			for (k = 0; k < j;) {
				if(k + BLOCK_SIZE < j){
					sum += A[IDX(i, k, n)] * A[IDX(j, k, n)] +
						A[IDX(i, k+1, n)] * A[IDX(j, k+1, n)] +
						A[IDX(i, k+2, n)] * A[IDX(j, k+2, n)] +
						A[IDX(i, k+3, n)] * A[IDX(j, k+3, n)] +
						A[IDX(i, k+4, n)] * A[IDX(j, k+4, n)] +
						A[IDX(i, k+5, n)] * A[IDX(j, k+5, n)] +
						A[IDX(i, k+6, n)] * A[IDX(j, k+6, n)] +
						A[IDX(i, k+7, n)] * A[IDX(j, k+7, n)];
					k += BLOCK_SIZE;
				}else{
					sum += A[IDX(i, k, n)] * A[IDX(j, k, n)];	
					k++;
				}
				
			}
			A[IDX(i, j, n)] -= sum;
			sum = 0.0;
		}

		if (A[IDX(j, j, n)] < 0.0) {
			return (1);
		}

		tmp = A[IDX(j, j, n)] = sqrt(A[IDX(j, j, n)]);
		for (i = j + 1; i < n;){
			if(i + BLOCK_SIZE < n){
				A[IDX(i, j, n)] /= tmp;
				A[IDX(i+1, j, n)] /= tmp;
				A[IDX(i+2, j, n)] /= tmp;
				A[IDX(i+3, j, n)] /= tmp;
				A[IDX(i+4, j, n)] /= tmp;
				A[IDX(i+5, j, n)] /= tmp;
				A[IDX(i+6, j, n)] /= tmp;
				A[IDX(i+7, j, n)] /= tmp;
				i += 8;
			}else{
				A[IDX(i, j, n)] /= A[IDX(j, j, n)];
				i++;
			}
		}
	}

	return (0);
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
	
	papi_init();
	
	A = generateSPDmatrix();
	
	papi_start();
	dtime = dclock();
 	chol_left_looking(A, SIZE);
 	dtime = dclock()-dtime;
 	papi_stop(0);
 	
	printf( "Time: %le \n", dtime);
	printf("papi\n");
	print_papi_results();

	return 0;
}
