#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "sparseMatrix.h"
#include "SpGEMV.h"
#include "parser.h"
#include "utils.h"
#include "macros.h"
#include "ompChunksDivide.h"
#include "ompGetICV.h"  //ICV - RUNTIME information audit auxs

#include <cblas.h>
double* SGEMVCBLAS(spmat* mat, double* inVect){
    CBLAS_LAYOUT layout=CblasRowMajor;
    CBLAS_TRANSPOSE notrans=CblasNoTrans;
    CBLAS_INT m=mat->M, n=mat->N;
    double* denseMat = CSRToDense(mat);
    if (!denseMat){
        ERRPRINT("GEMVCheckCBLAS: aux dense matrix alloc failed\n");
        return NULL;
    }
    double* oracleOut = malloc(m * sizeof(*oracleOut));
    if (!oracleOut){
        ERRPRINT("GEMVCheckCBLAS: out for serial oracle malloc failed\n");
        goto _err;
    }
    VERBOSE printf("computimg Sparse GEMV using densification over LAPACK.CBLAS\n");
    cblas_dgemv(layout,notrans,m,n,1.0,denseMat,n,inVect,1,0.0,oracleOut,1);
    
    _free:
    free(denseMat);
    return oracleOut;

    _err:
    if (denseMat)   free(denseMat);
    if (oracleOut)  free(oracleOut);
    return NULL;
}

