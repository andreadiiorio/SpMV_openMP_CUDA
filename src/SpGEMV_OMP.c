#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "SpGEMV.h"
#include "sparseMatrix.h"
#include "macros.h"
#include "utils.h"
#include "config.h"
#include "ompChunksDivide.h"

//global vars	->	audit
double Start,End,Elapsed,ElapsedInternal;

int spgemvRowsBasic(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;

    double acc;
    ulong c;
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc) (mat->M,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();

    #pragma omp parallel for schedule(runtime) private(acc)
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
    DEBUG printf("spgemvRowsBasic with %lu x %lu- %lu NZ CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,mat->NZ,ElapsedInternal);
    out = EXIT_SUCCESS;
    return out;
}

int spgemvRowsBlocks(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc;

    ulong rowBlock = mat->M / cfg->gridRows, rowBlockRem = mat->M % cfg->gridRows; 
    ulong block,startRow,c;
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc)(cfg->gridRows,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();

    #pragma omp parallel for schedule(runtime) private(acc, block, startRow)
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
    DEBUG printf("spgemvRowsBasic with %lu x %lu- %lu NZ CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,mat->NZ,ElapsedInternal);
    out = EXIT_SUCCESS;
    return out;
}

/*
 * spgemv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * 2D tiles of @mat are accessed via offsets from an aux allocated matrix map
 * final reduction to sum up all intermediate tiles' results
 */
int spgemvTiles(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc, *tilesOutTmp=NULL;
    ulong* offsets = NULL;   //aux matrix:   for each row -> colGroupStarts
    //2D INDEXING AUX
    ulong gridSize = cfg->gridRows * cfg->gridCols;
    ulong _rowBlock = mat->M / cfg->gridRows, _rowBlockRem = mat->M % cfg->gridRows; 
    //ulong _colBlock = mat->N / cfg->gridCols, _colBlockRem = mat->N % cfg->gridCols;
    ulong startRow,rowBlock; //,startCol,colBlock;
    if (!(offsets = colsOffsetsPartitioningUnifRanges(mat,cfg->gridCols))) goto _free;
    if (!(tilesOutTmp = malloc (mat->M * cfg->gridCols * sizeof(*tilesOutTmp)))){
        ERRPRINT("spgemvTiles:  tilesOutTmp malloc errd\n");
        goto _free;
    }
    memset(outVect,0,mat->M * sizeof(*outVect));
    ulong tileID,t_i,t_j,c;                            //for aux vars
    
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc)(gridSize,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();
    #pragma omp parallel for schedule(runtime) private(acc, rowBlock, startRow)
    for (tileID = 0; tileID < gridSize; tileID++){
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
    if (tilesOutTmp)    free(tilesOutTmp);
    return out;
}


/*
 * spgemv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * 2D tiles of @mat are accessed via allocated matrixes from column partitions of @mat
 * final reduction to sum up all intermediate tiles' results
 */
int spgemvTilesAllocd(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
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
        ERRPRINT("spgemvTiles:  tilesOutTmp malloc errd\n");
        goto _free;
    }
    memset(outVect,0,mat->M * sizeof(*outVect));
    
    ulong tileID,t_i,t_j,c;                            //for aux vars
    ((CHUNKS_DISTR_INTERF) cfg->chunkDistrbFunc)(gridSize,mat,cfg);
    AUDIT_INTERNAL_TIMES    Start = omp_get_wtime();
    #pragma omp parallel for schedule(runtime) private(acc, rowBlock, startRow)
    for (tileID = 0; tileID < gridSize; tileID++){
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
    for (ulong i=0; i<cfg->gridCols; i++)    freeSpmatInternal(colParts+i); //TODO FIRST ALLOCATION
    if (tilesOutTmp)    free(tilesOutTmp);
    return out;
}


int sgemvSerial(spmat* mat,double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc,end,start;
    start = omp_get_wtime();

    ulong i,j,k;
    omp_set_num_threads(cfg->threadNum);
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
