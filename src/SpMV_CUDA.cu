/*
 * 	dev:	Andrea Di Iorio
 * 	CUDA implementations of SpMV with CSR and ELL spMatrix rappresetantations
 */

extern "C" {
#include "sparseMatrix.h"
#include "SpGEMV.h"
}
#include "cudaUtils.h"

/*
 * fill @vShrd vector in shared mem with corrisponding values in @v 
 * using scheduled thread in a 1D grid of 1D blocks
 * 	-> use every thread available, reiterating if the whole grid is not enough
 */
__device__ void fillSharedVector_block1D_grid1D(double* v, double* vShrd, ulong len){ //TODO inline
	uint tid = threadIdx.x + blockIdx.x*blockDim.x;	//1D block,1D grid
	for(ulong i = tid; i < len; i += blockDim.x*gridDim.x)		vShrd[i] = v[i];
}
/*	template fill
	//copy v:	GMEM -> SHMEM because the lot of accesses it'll have
	extern __shared__ double* vShrd;	
	//TODO fit check (with concurr affect) + add vShrdVersion
	fillSharedVector_block1D_grid1D(v,vShrd,m->N);
	if (tRow < m->N)	vShrd[tRow] = v[tRow]; not horiz. rectangulare matrixes
	__syncthreads();
 */


__global__ void cudaSpMVRowsCSR(spmat* m,double* v,CONFIG cfg,double* outV){
	uint tRow = threadIdx.x + blockIdx.x*blockDim.x;	//1D block,1D grid
	if (tRow >= m->M)	return;
	ulong eIdx,rLen;
	#ifdef ROWLENS
	rLen = m->RL[tRow];
	#else
	rLen = m->IRP[tRow+1] - m->IRP[tRow];
	#endif
	double outEntry = 0;	//thread's output entry in outV
	for(ulong i=0; i<rLen; i++){
		eIdx = m->IRP[tRow] + i;
		outEntry += m->AS[eIdx] * v[m->JA[eIdx]];
	}
	outV[tRow] = outEntry;
}
__global__ void cudaSpMVWarpPerRowCSR(spmat* m,double* v,CONFIG cfg,double* outV){
}
__global__ void cudaSpMVRowsELL(spmat* m,double* v,CONFIG cfg,double* outV){
}
__global__ void cudaSpMVWarpPerRowELL(spmat* m,double* v,CONFIG cfg,double* outV){
}
__global__ void cudaSpMVWarpPerRowELLNTrasposed(spmat* m,double* v,CONFIG cfg,double* outV){
}
