#ifndef _SPGEMV
#define _SPGEMV

#include <assert.h>
#include "macros.h"
#include "sparseMatrix.h"

///COMPUTE MODES
#define CSR_SORTED_ROWS     "SORTED_ROWS"   ///TODO PARTITIONING ADD? 
#define CSR_ROWS            "CSR_ROWS"            
#define CSR_ROWS_GROUPS     "CSR_ROWS_GROUPS"
#define CSR_TILES           "CSR_TILES"           
#define CSR_TILES_ALLOCD    "CSR_TILES_ALLOCD"

#define ELL_ROWS			"ELL_ROWS"
#define ELL_ROWS_GROUPS		"ELL_ROWS_GROUPS"
#define ELL_TILES			"ELL_TILES"

#define CUDA_CSR_ROWS		"CUDA_CSR_ROWS"
#define CUDA_CSR_ROWS_WARP	"CUDA_CSR_ROWS_WARP"
#define CUDA_ELL_ROWS		"CUDA_ELL_ROWS"
#define CUDA_ELL_ROWS_WARP	"CUDA_ELL_ROWS_WARP"
#define CUDA_ELL_ROWS_WARP_NT	"CUDA_ELL_ROWS_WARP_NN_TRANSPOSED"
typedef enum {   
	//OMP MODES
	_CSR_SORTED_ROWS,
	_CSR_ROWS,
	_CSR_ROWS_GROUPS,
	_CSR_TILES,
	_CSR_TILES_ALLOCD,
	//ELL MODES
	_ELL_ROWS,
	_ELL_ROWS_GROUPS,
	_ELL_TILES,
	//CUDA MODES
	_CUDA_CSR_ROWS,
	_CUDA_CSR_ROWS_WARP,
	_CUDA_ELL_ROWS,
	_CUDA_ELL_ROWS_WARP,
	_CUDA_ELL_ROWS_WARP_NT,
} COMPUTE_MODE;  

#include "config.h"
//Sparse parallel Matrix-Vector computation @mat,@inVector,@config,@outVect
typedef int (SPGEMV)         		(spmat*,double*,CONFIG*,double*);
typedef int (*SPGEMV_INTERF) 		(spmat*,double*,CONFIG*,double*);

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

#ifdef __CUDACC__
typedef void (SPGEMV_CUDA)   		(spmat*,double*,CONFIG,double*);
typedef void (*SPGEMV_CUDA_INTERF)	(spmat*,double*,CONFIG,double*);

//__global__ void cudaSpMVRowsCSR(spmat* m,double* v,CONFIG cfg,double* outV);
__global__ SPGEMV_CUDA cudaSpMVRowsCSR;
__global__ SPGEMV_CUDA cudaSpMVWarpPerRowCSR;
__global__ SPGEMV_CUDA cudaSpMVRowsELL;
__global__ SPGEMV_CUDA cudaSpMVWarpPerRowELL;
__global__ SPGEMV_CUDA cudaSpMVWarpPerRowELLNTrasposed;

///threading sizeing
#define BLOCKS_1D	( 1u << 10 )
#endif
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
