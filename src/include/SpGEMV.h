#ifndef _SPGEMV
#define _SPGEMV

#include <assert.h>
#include "macros.h"
#include "sparseMatrix.h"

#define _ROWS            "ROWS"            
#define _SORTED_ROWS     "SORTED_ROWS"     
#define _TILES           "TILES"           
  
typedef enum {   
    ROWS,        
    SORTED_ROWS, 
    TILES        
} COMPUTE_MODE;  

#include "config.h"
//Sparse parallel Matrix-Vector computation @mat,@inVector,@config,@outVect
typedef int (SPGEMV)         (spmat*,double*,CONFIG*,double*);
typedef int (*SPGEMV_INTERF) (spmat*,double*,CONFIG*,double*);
/*
 * basic spgemv with row partitioning, 1 row per thread in consecutive order 
 */
//int spgemvRowsBasic(spmat* mat, double* vect, CONFIG* cfg, double* outVect);
SPGEMV spgemvRowsBasic;

/*
 * spgemv with row partitioning in blocks, consecutive order 
 */
SPGEMV spgemvRowsBlocks;

/*
 * spgemv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * 2D tiles of @mat are accessed via offsets from an aux allocated matrix map
 * final reduction to sum up all intermediate tiles' results
 */
SPGEMV spgemvTiles;
/*
 * spgemv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * 2D tiles of @mat are accessed via allocated matrixes from column partitions of @mat
 * final reduction to sum up all intermediate tiles' results
 */
SPGEMV spgemvTilesAllocd;

//serial implementaion of sparse GEMV for DEBUG
SPGEMV sgemvSerial;


static const SPGEMV_INTERF  SpgemvFuncs[] = {
    &sgemvSerial,
    &spgemvRowsBasic,
    &spgemvRowsBlocks,
    &spgemvTiles,
    &spgemvTilesAllocd,
    NULL,
};
#endif
