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

////////////////////////////////////// CSR /////////////////////////////////////
/*
 * basic spgemv with row partitioning, 1 row per thread in consecutive order 
 */
SPGEMV spgemvRowsBasicCSR;

/*
 * spgemv with rows partitioning in blocks, consecutive order 
 */
SPGEMV spgemvRowsBlocksCSR;

/*
 * spgemv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * 2D tiles of @mat are accessed via offsets from an aux allocated matrix map
 * final reduction to sum up all intermediate tiles' results
 */
SPGEMV spgemvTilesCSR;
/*
 * spgemv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * 2D tiles of @mat are accessed via allocated matrixes from column partitions of @mat
 * final reduction to sum up all intermediate tiles' results
 */
SPGEMV spgemvTilesAllocdCSR;
////////////////////////////////////// ELL /////////////////////////////////////
/* Basic ELL: MAX_ROW_NZ,AS,JA ==> no info to skip padding vals in omp (only cancel)
 * Defining ROWLENS -> mat->RL will be added 
 * -> easily determined exact end of each compute loop
 */

/*
 * basic spgemv with row partitioning, 1 row per thread in consecutive order 
 */
SPGEMV spgemvRowsBasicELL;

/*
 * spgemv with rows partitioning in blocks, consecutive order 
 */
SPGEMV spgemvRowsBlocksELL;

/*
 * spgemv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * final reduction to sum up all intermediate tiles' results
 */
SPGEMV spgemvTilesELL;
////////////////////////////////////////////////////////////////////////////////

//serial implementaion of sparse GEMV for DEBUG
SPGEMV sgemvSerial;

//wrap implementation func pointers for tests
static const SPGEMV_INTERF  SpgemvCSRFuncs[] = {
    &sgemvSerial,
    &spgemvRowsBasicCSR,
    &spgemvRowsBlocksCSR,
    &spgemvTilesCSR,
    &spgemvTilesAllocdCSR,
};
static const SPGEMV_INTERF  SpgemvELLFuncs[] = {
    &spgemvRowsBasicELL,
    &spgemvRowsBlocksELL,
    &spgemvTilesELL,
};
#endif
