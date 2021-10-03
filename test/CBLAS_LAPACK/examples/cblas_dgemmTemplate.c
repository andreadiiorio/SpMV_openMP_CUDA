/*
 * enahanced from https://gist.github.com/xianyi/5780018 ADAPTED FOR LAPACK.CBLAS to perform a regular matrix multiplication
 *
 * BUILDABLE FROM CBLAS DIR -- STATIC LINK WITH gfortran
 *  gcc -O3 -I../include -c -o cblas_dgemmTemplate.o cblas_dgemmTemplate.c
 *  gfortran -O2 -frecursive  -o cblas_dgemmTemplate  cblas_dgemmTemplate.o  ../../libcblas.a ../../librefblas.a
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "cblas.h"

typedef unsigned int uint;
static void printMatrix(uint M,uint N,double mat[][N]);
int main()
{
    CBLAS_LAYOUT Layout=CblasRowMajor;;
    CBLAS_TRANSPOSE transa=CblasNoTrans, transb=CblasNoTrans;

    uint i,m=4, n=6, k=5, lda=k, ldb=n, ldc=n;

	uint sizeofa = m * k;
	uint sizeofb = k * n;
	uint sizeofc = m * n;
	double alpha = 1, beta = 1;


	double* A = malloc(sizeof(*A) * sizeofa);
	double* B = malloc(sizeof(*B) * sizeofb);
	double* C = malloc(sizeof(*C) * sizeofc);

	//srand((unsigned)time(NULL));
	for (i=0; i<sizeofa; i++)		A[i] = i%3+1;//(rand()%100)/10.0;
	for (i=0; i<sizeofb; i++)		B[i] = i%3+1;//(rand()%100)/10.0;
	//for (i=0; i<sizeofc; i++)		C[i] = i%3+1;//(rand()%100)/10.0;
	printMatrix(m,k,(double (*)[k]) A); printMatrix(k,n,(double (*)[n])B);
	printf("dgemm:: m=%d,n=%d,k=%d,alpha=%lf,beta=%lf,sizeofc=%d\n",m,n,k,alpha,beta,sizeofc);
	
    struct timeval start,finish;double duration;
	gettimeofday(&start, NULL);
	////// CBLAS CALL
    cblas_dgemm( Layout, transa,transb, m, n, k, alpha, A, lda, B, ldb, beta,C, ldc );
	///////////////////
	gettimeofday(&finish, NULL);
    printMatrix(m,n,(double (*)[n])C);

    free(A);free(B);free(C);
    return 0;
}
static void printMatrix(uint M,uint N,double mat[][N]){
    for(uint i=0; i<M; i++){
        for(uint j=0; j<N; j++)     printf("%1.1le ",mat[i][j]);
        printf("\n");
    }
    printf("\n\n");
}

