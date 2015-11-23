#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "emmintrin.h"

/* for Fortran - j*n + i */
//#define IDX(i, j, n)	(((i) * (n)) + (j))
#define IDX(i, j, n) (((j)+ (i)*(n)))
#define BLOCK_SIZE 16
#define SIZE 16

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
			printf("%le \t", tmp);
			A[IDX(i, j, SIZE)] = tmp;
			A[IDX(j, i, SIZE)] = tmp;
		}
	}

	mm(A, At, SPDM);

	return SPDM;
}

int chol_left_looking(double *A, unsigned int n){
	register unsigned int i, j, k;
	register double sum, tmp;
	register __m128d tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8,
					 tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15;
	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = A[IDX(i, j, n)];
			for (k = 0; k < j; ) {
				if(k + BLOCK_SIZE < j){
					//load data
					tmp0 = _mm_loadu_pd(A+IDX(i, k, n)); 
					tmp1 = _mm_loadu_pd(A+IDX(j, k, n));
					tmp2 = _mm_loadu_pd(A+IDX(i, k+2, n));
					tmp3 = _mm_loadu_pd(A+IDX(j, k+2, n));
					tmp4 = _mm_loadu_pd(A+IDX(i, k+4, n));
					tmp5 = _mm_loadu_pd(A+IDX(j, k+4, n));
					tmp6 = _mm_loadu_pd(A+IDX(i, k+6, n));
					tmp7 = _mm_loadu_pd(A+IDX(j, k+6, n));
					tmp8 = _mm_loadu_pd(A+IDX(i, k+8, n));
					tmp9 = _mm_loadu_pd(A+IDX(j, k+8, n));
					tmp10 = _mm_loadu_pd(A+IDX(i, k+10, n));
					tmp11 = _mm_loadu_pd(A+IDX(j, k+10, n));
					tmp12 = _mm_loadu_pd(A+IDX(i, k+12, n));
					tmp13 = _mm_loadu_pd(A+IDX(j, k+12, n));
					tmp14 = _mm_loadu_pd(A+IDX(i, k+14, n));
					tmp15 = _mm_loadu_pd(A+IDX(j, k+14, n));

					//multiplication A[IDX(i, k, n)] * A[IDX(j, k, n)] and A[IDX(i, k+1, n)] * A[IDX(j, k+1, n)]
					tmp0 = _mm_mul_pd(tmp0, tmp1); 
					tmp2 = _mm_mul_pd(tmp2, tmp3);
					tmp4 = _mm_mul_pd(tmp4, tmp5);
					tmp6 = _mm_mul_pd(tmp6, tmp7);
					tmp8 = _mm_mul_pd(tmp8, tmp9);
					tmp10 = _mm_mul_pd(tmp10, tmp11);
					tmp12 = _mm_mul_pd(tmp12, tmp13);
					tmp14 = _mm_mul_pd(tmp14, tmp15);


					/*
					sum value like binary tree 
													  tmp0
									  /                               \
									tmp0                             tmp8
							  /	             \                /               \
							tmp0		    tmp4            tmp8            tmp12 
 						  /		 \        /	     \        /      \        /      \
						tmp0	tmp2	tmp4	tmp6	tmp8	tmp10	tmp12	tmp14				
					*/
					tmp0 = _mm_add_pd(tmp0, tmp2); 
					tmp4 = _mm_add_pd(tmp4, tmp6);
					tmp8 = _mm_add_pd(tmp8, tmp10);
					tmp12 = _mm_add_pd(tmp12, tmp14);

					tmp0 = _mm_add_pd(tmp0, tmp4);
					tmp8 = _mm_add_pd(tmp8, tmp12);

					tmp0 = _mm_add_pd(tmp0, tmp8);

					sum -= tmp0[0] + tmp0[1];
					
					k += BLOCK_SIZE;
				}else if(k + (BLOCK_SIZE / 2) < j){
					tmp0 = _mm_loadu_pd(A+IDX(i, k, n)); 
					tmp1 = _mm_loadu_pd(A+IDX(j, k, n));
					tmp2 = _mm_loadu_pd(A+IDX(i, k+2, n));
					tmp3 = _mm_loadu_pd(A+IDX(j, k+2, n));
					tmp4 = _mm_loadu_pd(A+IDX(i, k+4, n));
					tmp5 = _mm_loadu_pd(A+IDX(j, k+4, n));
					tmp6 = _mm_loadu_pd(A+IDX(i, k+6, n));
					tmp7 = _mm_loadu_pd(A+IDX(j, k+6, n));

					tmp0 = _mm_mul_pd(tmp0, tmp1); 
					tmp2 = _mm_mul_pd(tmp2, tmp3);
					tmp4 = _mm_mul_pd(tmp4, tmp5);
					tmp6 = _mm_mul_pd(tmp6, tmp7);

					tmp0 = _mm_add_pd(tmp0, tmp2); 
					tmp4 = _mm_add_pd(tmp4, tmp6);

					tmp0 = _mm_add_pd(tmp0, tmp4);

					sum -= tmp0[0] + tmp0[1];

					k += BLOCK_SIZE / 2;

				}else if(k + (BLOCK_SIZE / 4) < j){
					tmp0 = _mm_loadu_pd(A+IDX(i, k, n)); 
					tmp1 = _mm_loadu_pd(A+IDX(j, k, n));
					tmp2 = _mm_loadu_pd(A+IDX(i, k+2, n));
					tmp3 = _mm_loadu_pd(A+IDX(j, k+2, n));

					tmp0 = _mm_mul_pd(tmp0, tmp1); 
					tmp2 = _mm_mul_pd(tmp2, tmp3);
					
					tmp0 = _mm_add_pd(tmp0, tmp2); 

					sum -= tmp0[0] + tmp0[1];

					k += BLOCK_SIZE / 4;
				}else if(k + (BLOCK_SIZE / 8) < j){
					tmp0 = _mm_loadu_pd(A+IDX(i, k, n)); 
					tmp1 = _mm_loadu_pd(A+IDX(j, k, n));

					tmp0 = _mm_mul_pd(tmp0, tmp1); 

					sum -= tmp0[0] + tmp0[1];

					k += BLOCK_SIZE / 8;
				}else{ 
					sum -= A[IDX(i, k, n)] * A[IDX(j, k, n)];	
					k++;
				}
				
			}
			A[IDX(i, j, n)] = sum;
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

	return 0;
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
	int i, j;
	A = generateSPDmatrix();
	for(i = 0; i < SIZE; i++){
		for(j = 0; j < SIZE; j++){
			printf("%le \t", A[IDX(i, j, SIZE)]);	
		}
		printf("\n");
	}
	dtime = dclock();
 	chol_left_looking(A, SIZE);
 	dtime = dclock()-dtime;
 	double gflops = ((1.0/3.0) * SIZE * SIZE * SIZE * 10e-9) / dtime;
	printf( "Time: %le \n", dtime);
	printf("Gflops: %le \n", gflops);

	

	
	return 0;
}
