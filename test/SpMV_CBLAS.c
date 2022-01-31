/*
Copyright Andrea Di Iorio 2022
This file is part of SpMV_OMP_CUDA
SpMV_OMP_CUDA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SpMV_OMP_CUDA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SpMV_OMP_CUDA.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "sparseMatrix.h"
#include "SpMV.h"
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

