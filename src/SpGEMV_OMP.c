#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "SpGEMV.h"
#include "sparseMatrix.h"
#include "macros.h"
#include "utils.h"

//global vars	->	audit
double Start,End,Elapsed,ElapsedInternal;

int spgemvRowsBasic(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc;
    Start = omp_get_wtime();

    uint i,j,k;
    #pragma omp parallel for schedule(static) shared(vect,mat) private(i, j, k)
    for (i=0;i<mat->M;i++){
        acc = 0;
        for (j=mat->IRP[i]; j<mat->IRP[i+1]; j++){
           k = mat->JA[j];
           acc += mat->AS[j] * vect[k];
        }
        outVect[i] = acc;
    }
 
    End = omp_get_wtime();
    Elapsed = ElapsedInternal = End - Start;
    VERBOSE printf("spgemvRowsBasic with %u x %u- %u NZ CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,mat->NZ,ElapsedInternal);
    out = EXIT_SUCCESS;
    return out;
}

int sgemvSerial(spmat* mat,double* vect, CONFIG* cfg, double* outVect){
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
    VERBOSE printf("sgemvSerial with %u x %u - %u NNZ  CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,mat->NZ,end-start);
    out = EXIT_SUCCESS;
    return out;
}

