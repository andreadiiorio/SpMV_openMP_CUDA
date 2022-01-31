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
//ELL sparse matrixes SpMV implementation for OMP

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

int spmvRowsBasicELL(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;

    double acc;
    ulong c,rMax = mat->MAX_ROW_NZ,rowStart,rowEnd;
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc) (mat->M,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();

    #pragma omp parallel for schedule(runtime) private(acc,c,rowStart,rowEnd)
    for (ulong r=0;  r<mat->M;   r++){
        acc = 0;
        
        rowStart    = IDX2D(r,0,rMax);
        rowEnd      = IDX2D(r,rMax,rMax);
        #ifdef ROWLENS
        rowEnd      = rowStart+mat->RL[r];
        #endif
        #if SIMD_ROWS_REDUCTION == TRUE
        #pragma omp simd reduction(+:acc)
        #endif
        for (ulong j = rowStart; j < rowEnd; j++){
            c = mat->JA[j];
            acc += mat->AS[j] * vect[c];
        }

        outVect[r] = acc;
    }
 
    AUDIT_INTERNAL_TIMES{
        End=omp_get_wtime();
        ElapsedInternal = End-Start;
    }
    out = EXIT_SUCCESS;
    return out;
}

int spmvRowsBlocksELL(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc;

    ulong rowBlock = mat->M / cfg->gridRows, rowBlockRem = mat->M % cfg->gridRows; 
    ulong block,startRow,c,rMax = mat->MAX_ROW_NZ,rowStart,rowEnd;
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc)(cfg->gridRows,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();

    #pragma omp parallel for schedule(runtime) private(acc,block,startRow,c,rowStart,rowEnd)
    for (ulong b=0;   b<cfg->gridRows;   b++){
        block      = UNIF_REMINDER_DISTRI(b,rowBlock,rowBlockRem);
        startRow   = UNIF_REMINDER_DISTRI_STARTIDX(b,rowBlock,rowBlockRem);
        for (ulong r=startRow;  r<startRow+block;  r++){
            acc = 0;
        
            rowStart    = IDX2D(r,0,rMax);
            rowEnd      = IDX2D(r,rMax,rMax);
            #ifdef ROWLENS
            rowEnd      = rowStart+mat->RL[r];
            #endif
            #if SIMD_ROWS_REDUCTION == TRUE
            #pragma omp simd reduction(+:acc)
            #endif
            for (ulong j = rowStart; j < rowEnd; j++){
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
    out = EXIT_SUCCESS;
    return out;
}

int spmvTilesELL(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc, *tilesOutTmp=NULL;
    //2D INDEXING AUX
    ulong gridSize = cfg->gridRows * cfg->gridCols,rMax = mat->MAX_ROW_NZ;
    ulong _rowBlock = mat->M / cfg->gridRows, _rowBlockRem = mat->M % cfg->gridRows; 
    ulong _colBlock = rMax / cfg->gridCols, _colBlockRem = rMax % cfg->gridCols;
    ulong startRow,rowBlock,startCol,colBlock,  rowPartStart,rowPartEnd;
    if (!(tilesOutTmp = malloc(mat->M * cfg->gridCols * sizeof(*tilesOutTmp)))){
        ERRPRINT("spmvTiles:  tilesOutTmp malloc errd\n");
        goto _free;
    }
    memset(outVect,0,mat->M * sizeof(*outVect));
    ulong t_i,t_j,c;                            //for aux vars
    
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc)(gridSize,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();
    #pragma omp parallel for schedule(runtime) \
      private(acc,rowBlock,startRow,startCol,colBlock,c,t_i,t_j,rowPartStart,rowPartEnd)
    for (ulong tileID = 0; tileID < gridSize; tileID++){
        ///get iteration's indexing variables
        //tile index in the 2D grid of AB computation TODO OMP HOW TO PARALLELIZE 2 FOR
        t_i = tileID/cfg->gridCols;  //i-th row block
        t_j = tileID%cfg->gridCols;  //j-th col block
        //get tile row-cols group FAIR sizes
        rowBlock = UNIF_REMINDER_DISTRI(t_i,_rowBlock,_rowBlockRem); 
        startRow = UNIF_REMINDER_DISTRI_STARTIDX(t_i,_rowBlock,_rowBlockRem);
        colBlock = UNIF_REMINDER_DISTRI(t_j,_colBlock,_colBlockRem); 
        startCol = UNIF_REMINDER_DISTRI_STARTIDX(t_j,_colBlock,_colBlockRem);

        for (ulong r=startRow; r<startRow+rowBlock;  r++){
            acc = 0;
            rowPartStart = IDX2D(r,startCol,rMax);
            rowPartEnd   = rowPartStart + colBlock; 
            #ifdef ROWLENS
            rowPartEnd   = MIN(rowPartEnd,IDX2D(r,0,rMax) + mat->RL[r]);
            #endif
            #if SIMD_ROWS_REDUCTION == TRUE
            #pragma omp simd reduction(+:acc)
            #endif
            //2D CSR tile SpMV using cols partition of row r
            for (ulong j=rowPartStart; j<rowPartEnd; j++){
                c = mat->JA[j];
                acc += mat->AS[j] * vect[c]; 
            }
            tilesOutTmp[ IDX2D(r,t_j,cfg->gridCols) ] = acc;
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
    if(tilesOutTmp)     free(tilesOutTmp);
    return out;
}

