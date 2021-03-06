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

//dev:  Andrea Di Iorio - 0277550
//CSR sparse matrixes SpMV implementation for OMP

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "SpMV.h"
#include "sparseMatrix.h"
#include "macros.h"
#include "utils.h"
#include "config.h"
#include "ompChunksDivide.h"


int spmvRowsBasicCSR(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;

    double acc;
    ulong c;
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc) (mat->M,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();

    #pragma omp parallel for schedule(runtime) private(acc,c)
    for (ulong r=0;  r<mat->M;   r++){
        acc = 0;
        #if SIMD_ROWS_REDUCTION == TRUE
        #pragma omp simd reduction(+:acc)
        #endif
        for (ulong j=mat->IRP[r]; j<mat->IRP[r+1]; j++){
           c = mat->JA[j];
           acc += mat->AS[j] * vect[c];
        }
        outVect[r] = acc;
    }
 
    AUDIT_INTERNAL_TIMES{
        End=omp_get_wtime();
        ElapsedInternal = End-Start;
    }
    DEBUG printf("spmvRowsBasic with %lu x %lu- %lu NZ CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,mat->NZ,ElapsedInternal);
    out = EXIT_SUCCESS;
    return out;
}

int spmvRowsBlocksCSR(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc;

    ulong rowBlock = mat->M / cfg->gridRows, rowBlockRem = mat->M % cfg->gridRows; 
    ulong block,startRow,c;
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc)(cfg->gridRows,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();

    #pragma omp parallel for schedule(runtime) private(acc,block,startRow,c)
    for (ulong b=0;   b<cfg->gridRows;   b++){
        block      = UNIF_REMINDER_DISTRI(b,rowBlock,rowBlockRem);
        startRow   = UNIF_REMINDER_DISTRI_STARTIDX(b,rowBlock,rowBlockRem);
        for (ulong r=startRow;  r<startRow+block;  r++){
            acc = 0;
            #if SIMD_ROWS_REDUCTION == TRUE
            #pragma omp simd reduction(+:acc)
            #endif
            for (ulong j=mat->IRP[r]; j<mat->IRP[r+1]; j++){
               c = mat->JA[j];
               acc += mat->AS[j] * vect[c];
            }
            outVect[r] = acc;
        }
    }
    AUDIT_INTERNAL_TIMES{
        End=omp_get_wtime();
        ElapsedInternal = End-Start;
    }
    DEBUG printf("spmvRowsBasic with %lu x %lu- %lu NZ CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,mat->NZ,ElapsedInternal);
    out = EXIT_SUCCESS;
    return out;
}

int spmvTilesCSR(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc, *tilesOutTmp=NULL;
    ulong* offsets = NULL;   //aux matrix:   for each row -> colGroupStarts
    //2D INDEXING AUX
    ulong gridSize = cfg->gridRows * cfg->gridCols;
    ulong _rowBlock = mat->M / cfg->gridRows, _rowBlockRem = mat->M % cfg->gridRows; 
    //ulong _colBlock = mat->N / cfg->gridCols, _colBlockRem = mat->N % cfg->gridCols;
    ulong startRow,rowBlock; //,startCol,colBlock;
    if (!(offsets = colsOffsetsPartitioningUnifRanges(mat,cfg->gridCols))) goto _free;
    if (!(tilesOutTmp = malloc(mat->M * cfg->gridCols * sizeof(*tilesOutTmp)))){
        ERRPRINT("spmvTiles:  tilesOutTmp malloc errd\n");
        goto _free;
    }
    memset(outVect,0,mat->M * sizeof(*outVect));
    ulong t_i,t_j,c;                            //for aux vars
    
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc)(gridSize,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();
    #pragma omp parallel for schedule(runtime) private(acc,rowBlock,startRow,c,t_i,t_j)
    for (ulong tileID = 0; tileID < gridSize; tileID++){
        ///get iteration's indexing variables
        //tile index in the 2D grid of AB computation TODO OMP HOW TO PARALLELIZE 2 FOR
        t_i = tileID/cfg->gridCols;  //i-th row block
        t_j = tileID%cfg->gridCols;  //j-th col block
        //get tile row-cols group FAIR sizes
        rowBlock = UNIF_REMINDER_DISTRI(t_i,_rowBlock,_rowBlockRem); 
        startRow = UNIF_REMINDER_DISTRI_STARTIDX(t_i,_rowBlock,_rowBlockRem);
        //startCol = UNIF_REMINDER_DISTRI_STARTIDX(t_j,_colBlock,_colBlockRem);

        for (ulong r=startRow,partOffID;  r<startRow+rowBlock;  r++){
            partOffID = IDX2D(r,t_j,cfg->gridCols);
            acc = 0;
            #if SIMD_ROWS_REDUCTION == TRUE
            #pragma omp simd reduction(+:acc)
            #endif
            //2D CSR tile SpMV using cols partition of row r
            for (ulong j=offsets[partOffID]; j<offsets[partOffID+1]; j++){
                c = mat->JA[j];
                acc += mat->AS[j] * vect[c]; 
            }
            tilesOutTmp[partOffID] = acc;
        }
    }
    //#pragma omp parallel for reduction(+:outVect[:mat->M])    //TODO RARE SEGFAULT
    for (ulong r=0;  r<mat->M;   r++){
        for (ulong p=0;   p<cfg->gridCols;   p++){
            outVect[r] += tilesOutTmp[ IDX2D(r,p,cfg->gridCols) ];
        }
    }
    AUDIT_INTERNAL_TIMES{
        End=omp_get_wtime();
        ElapsedInternal = End-Start;
    }
    
    out = EXIT_SUCCESS;
    _free:
    if(offsets)         free(offsets);
    if(tilesOutTmp)     free(tilesOutTmp);
    return out;
}


int spmvTilesAllocdCSR(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc, *tilesOutTmp=NULL;
    spmat *colParts = NULL, *colPart;
    //2D INDEXING AUX
    ulong gridSize = cfg->gridRows * cfg->gridCols;
    ulong _rowBlock = mat->M / cfg->gridRows, _rowBlockRem = mat->M % cfg->gridRows; 
    //ulong _colBlock = mat->N / cfg->gridCols, _colBlockRem = mat->N % cfg->gridCols;
    ulong startRow,rowBlock; //,startCol,colBlock;
    if (!(colParts = colsPartitioningUnifRanges(mat,cfg->gridCols)))   return EXIT_FAILURE;
    if (!(tilesOutTmp = malloc (mat->M * cfg->gridCols * sizeof(*tilesOutTmp)))){
        ERRPRINT("spmvTiles:  tilesOutTmp malloc errd\n");
        goto _free;
    }
    memset(outVect,0,mat->M * sizeof(*outVect));
    
    ulong t_i,t_j,c;                            //for aux vars
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc)(gridSize,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();
    #pragma omp parallel for schedule(runtime) private(acc,rowBlock,startRow,c,t_i,t_j)
    for (ulong tileID = 0; tileID < gridSize; tileID++){
        ///get iteration's indexing variables
        //tile index in the 2D grid of AB computation TODO OMP HOW TO PARALLELIZE 2 FOR
        t_i = tileID/cfg->gridCols;  //i-th row block
        t_j = tileID%cfg->gridCols;  //j-th col block
        colPart = colParts + t_j;
        //get tile row-cols group FAIR sizes
        rowBlock = UNIF_REMINDER_DISTRI(t_i,_rowBlock,_rowBlockRem); 
        startRow = UNIF_REMINDER_DISTRI_STARTIDX(t_i,_rowBlock,_rowBlockRem);
        //startCol = UNIF_REMINDER_DISTRI_STARTIDX(t_j,_colBlock,_colBlockRem);
        for (ulong r=startRow,partOffID;  r<startRow+rowBlock;  r++){
            acc = 0;
            partOffID = IDX2D(r,t_j,cfg->gridCols);
            #if SIMD_ROWS_REDUCTION == TRUE
            #pragma omp simd reduction(+:acc)
            #endif
            for (ulong j=colPart->IRP[r]; j<colPart->IRP[r+1]; j++){
                c = colPart->JA[j];
                acc += colPart->AS[j] * vect[c]; 
            }
            tilesOutTmp[partOffID] = acc;
        }
    }
    //#pragma omp parallel for reduction(+:outVect[:mat->M])    //TODO RARE SEGFAULT
    for (ulong r=0;  r<mat->M;   r++){
        for (ulong p=0;   p<cfg->gridCols;   p++){
            outVect[r] += tilesOutTmp[ IDX2D(r,p,cfg->gridCols) ];
        }
    }
    AUDIT_INTERNAL_TIMES{
        End=omp_get_wtime();
        ElapsedInternal = End-Start;
    }
    
    out = EXIT_SUCCESS;
    _free:
    for (ulong i=0; i<cfg->gridCols; i++)    freeSpmatInternal(colParts+i); 
    free(colParts);
    if (tilesOutTmp)    free(tilesOutTmp);
    return out;
}


int sgemvSerial(spmat* mat,double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc,end,start;
    start = omp_get_wtime();

    ulong i,j,k;
    for (i=0;i<mat->M;i++){
        acc = 0;
        for (j=mat->IRP[i]; j<mat->IRP[i+1]; j++){
           k = mat->JA[j];
           acc += mat->AS[j] * vect[k];
        }
        outVect[i] = acc;
    }
 
    end = omp_get_wtime();
    DEBUG printf("sgemvSerial with %lu x %lu - %lu NNZ  CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,mat->NZ,end-start);
    out = EXIT_SUCCESS;
    return out;
}
