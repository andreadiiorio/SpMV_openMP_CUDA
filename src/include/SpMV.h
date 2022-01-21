#ifndef _SPMV
#define _SPMV

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
typedef int (SPMV)         		(spmat*,double*,CONFIG*,double*);
typedef int (*SPMV_INTERF) 		(spmat*,double*,CONFIG*,double*);

////////////////////////////////////// CSR /////////////////////////////////////
/*
 * basic spmv with row partitioning, 1 row per thread in consecutive order 
 */
SPMV spmvRowsBasicCSR;

/*
 * spmv with rows partitioning in blocks, consecutive order 
 */
SPMV spmvRowsBlocksCSR;

/*
 * spmv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * 2D tiles of @mat are accessed via offsets from an aux allocated matrix map
 * final reduction to sum up all intermediate tiles' results
 */
SPMV spmvTilesCSR;
/*
 * spmv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * 2D tiles of @mat are accessed via allocated matrixes from column partitions of @mat
 * final reduction to sum up all intermediate tiles' results
 */
SPMV spmvTilesAllocdCSR;
////////////////////////////////////// ELL /////////////////////////////////////
/* Basic ELL: MAX_ROW_NZ,AS,JA ==> no info to skip padding vals in omp (only cancel)
 * Defining ROWLENS -> mat->RL will be added 
 * -> easily determined exact end of each compute loop
 */

/*
 * basic spmv with row partitioning, 1 row per thread in consecutive order 
 */
SPMV spmvRowsBasicELL;

/*
 * spmv with rows partitioning in blocks, consecutive order 
 */
SPMV spmvRowsBlocksELL;

/*
 * spmv via 2D decomposition of matrix @mat with @cfg
 * each tile of the matrix produce a part of the result in a tmp allocated space
 * final reduction to sum up all intermediate tiles' results
 */
SPMV spmvTilesELL;
////////////////////////////////////////////////////////////////////////////////

//serial implementaion of sparse GEMV for DEBUG
SPMV sgemvSerial;

#ifdef __CUDACC__
typedef void (SPMV_CUDA)   		(spmat*,double*,CONFIG,double*);
typedef void (*SPMV_CUDA_INTERF)	(spmat*,double*,CONFIG,double*);

//__global__ void cudaSpMVRowsCSR(spmat* m,double* v,CONFIG cfg,double* outV);
__global__ SPMV_CUDA cudaSpMVRowsCSR;
__global__ SPMV_CUDA cudaSpMVWarpPerRowCSR;
__global__ SPMV_CUDA cudaSpMVRowsELL;
__global__ SPMV_CUDA cudaSpMVRowsELLNNTransposed;
//__global__ SPMV_CUDA cudaSpMVWarpPerRowELL;
__global__ SPMV_CUDA cudaSpMVWarpsPerRowELLNTrasposed;
///func pntrs for auto. tests-perfG.
static const SPMV_CUDA_INTERF SpmvCUDA_CSRFuncs[] = {
	&cudaSpMVRowsCSR,
	&cudaSpMVWarpPerRowCSR  //SpmvCUDA_CSRFuncs_WarpPerRowIdx 
};
#define SpmvCUDA_CSRFuncs_WarpPerRowIdx 1

static const SPMV_CUDA_INTERF SpmvCUDA_ELLFuncs[] = {
	&cudaSpMVRowsELL,
	&cudaSpMVRowsELLNNTransposed,
	&cudaSpMVWarpsPerRowELLNTrasposed
};
#define SpmvCUDA_ELLFuncs_NN_TraposedImpl 	1 
#define SpmvCUDA_ELLFuncs_WarpPerRowIdx 	2

///threading sizeing
#define BLOCKS_1D			( 1u << 8 )
#define BLOCKS_2D_WARP_R	( 1u << (10-5) )
#endif //__CUDACC__


//wrap implementation func pointers for tests
static const SPMV_INTERF  SpmvCSRFuncs[] = {
    &sgemvSerial,
    &spmvRowsBasicCSR,
    &spmvRowsBlocksCSR,
    &spmvTilesCSR,
    &spmvTilesAllocdCSR,
};
static const SPMV_INTERF  SpmvELLFuncs[] = {
    &spmvRowsBasicELL,
    &spmvRowsBlocksELL,
    &spmvTilesELL,
};

#endif
