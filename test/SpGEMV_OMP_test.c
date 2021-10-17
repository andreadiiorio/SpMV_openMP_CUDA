#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include <cblas.h>

#include "SpGEMV.h"
#include "parser.h"
#include "utils.h"
#include <macros.h>
#include "sparseMatrix.h"

//TODO MOVE TO HEDER
int serialDenseGEMVTest(spmat* mat, double* inVect, double* outToCheck);
typedef struct{
    ushort gridRows;
    ushort gridCols;
    //TODO FULL CONFIG DOCCED HERE
    int threadNum;  //thread num to use in an OMP parallel region ...
} CONFIG;
CONFIG Conf = {
    .gridRows = 8,
};

//compute function interface and its pointer definitions
typedef int (COMPUTEFUNC) (spmat*,double*,CONFIG*,double*);
typedef int (*COMPUTEFUNC_INTERF) (spmat*,double*,CONFIG*,double*);
COMPUTEFUNC spgemvRowsBasic;

#define HELP "usage: MatrixMarket_sparse_matrix_COO, vectorFile || "RNDVECT"\n"
int main(int argc, char** argv){
    int out=EXIT_FAILURE;
    if (init_urndfd())  return out;
    if (argc < 3 )  {ERRPRINT(HELP); return out;}
    
    double *vector = NULL, *outVector = NULL;
    uint vectSize;
    spmat* mat = NULL; 
    ////parse sparse matrix and dense vector
    if (!(mat = MMtoCSR(argv[1]))){
        ERRPRINT("err during conversion MM -> CSR\n");
        return out;
    }
    ////get the vector
    if (!(strncmp(argv[2],RNDVECT,strlen(RNDVECT)))){ //generate a random vector
        vectSize = mat->N;  //size for GEMV
        if (!(vector = malloc(vectSize * sizeof(*vector)))){
            ERRPRINT("rnd vector malloc failed\n");
            goto _free;
        }
        if (fillRndVector(vectSize,vector)){
            ERRPRINT("fillRndVector errd\n");
            goto _free;
        }
    } else{ //read vector from the given file
        if (!(vector = readVector(argv[2],&vectSize))){
            fprintf(stderr,"err during readVector at:%s\n",argv[2]);
            goto _free;
        }
        CONSISTENCY_CHECKS{
            if (vectSize != mat->N){
                ERRPRINT("vector not compatible with sparse matrix\n");
                goto _free;
            }
        }
    }
    //alloc space for 2 output vectors
    if (!(outVector = malloc( 2 * mat->M * sizeof(*outVector)))){
        ERRPRINT("outVector malloc errd\n");
        goto _free;
    }
    
    DEBUG {
        printf("sparse matrix:\n");printSparseMatrix(mat,TRUE);
        printf("vector:\n");printVector(vector,vectSize);
    }
  
    //// PARALLEL COMPUTATIONs
    //int maxThreads = omp_get_max_threads();
    //all parallel implementations
    if ((out = spgemvRowsBasic(mat,vector,&Conf,outVector))){
        ERRPRINT("compute function failed...\n"); goto _free;
    }
    printf("outVector:\n");printVector(outVector,vectSize);
    //// SERIAL - LAPACK.CBLAS COMPUTATION
    if (serialDenseGEMVTest(mat,vector,outVector)){
        ERRPRINT("LAPACK.CBLAS SERIAL, DENSE GEMV TEST FAILED!!\n");
        goto _free;
    }
    _free:
    if (mat)          free(mat);
    if (vector)       free(vector);
    if (outVector)    free(outVector);
    return out;
}
int serialDenseGEMVTest(spmat* mat, double* inVect, double* outToCheck){
    int out = EXIT_FAILURE;
    CBLAS_LAYOUT layout=CblasRowMajor;
    CBLAS_TRANSPOSE notrans=CblasNoTrans;
    CBLAS_INT m=mat->M, n=mat->N;
    double* denseMat = CSRToDense(mat);
    if (!denseMat){
        ERRPRINT("serialDenseGEMVTest aux dense matrix alloc failed\n");
        goto _free;
    }
    double* oracleOut = malloc(m * sizeof(*oracleOut));
    if (!oracleOut){
        ERRPRINT("serialDenseGEMVTest out for serial oracle malloc failed\n");
        goto _free;
    }
    VERBOSE printf("checking parallel implementation using LAPACK.CBLAS\n");
    cblas_dgemv(layout,notrans,m,n,1.0,denseMat,n,inVect,1,0.0,oracleOut,1);
    out = doubleVectorsDiff(outToCheck,oracleOut,m);
    
    _free:
    free(denseMat);
    return out;
}
//TODO CATALOG COMPUTE FUNC IN SEP. FILES IF NEEDED
/*
 * basic spgemv with row partitioning, 1 row per thread in consecutive order
 * statically assigned to threads
 */
int spgemvRowsBasic(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc,end,start;
    start = omp_get_wtime();

    uint i,j,k;
    omp_set_num_threads(cfg->threadNum);
    #pragma omp parallel for schedule(static) shared(vect,mat) private(i, j, k)
    for (i=0;i<mat->M;i++){
        acc = 0;
        for (j=mat->IRP[i]; j<mat->IRP[i+1]; j++){
           k = mat->JA[j];
           acc += mat->AS[j] * vect[k];
        }
        outVect[i] = acc;
    }
 
    end = omp_get_wtime();
    VERBOSE printf("spgemvRowsBasic with %u x %u CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,end-start);
    out = EXIT_SUCCESS;
    return out;
}
