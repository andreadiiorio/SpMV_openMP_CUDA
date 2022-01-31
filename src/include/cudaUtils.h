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

#ifndef UTILSCUDA_H
#define UTILSCUDA_H

#define WARPSIZE	32	//dimensioning in non device code -> TODO  normal constant not available
#include "sparseMatrix.h"

////ERROR MNGM WRAP
///from cuda samples helper header...
template <typename T>
void check(T result, char const *const func, const char *const file,
           int const line) {
  if (result) {
    fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
            static_cast<unsigned int>(result), cudaGetErrorName(result), func);
    exit(EXIT_FAILURE);
  }
}
// This will output the proper CUDA error strings in the event
// that a CUDA host call returns an error
#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)
///
/*
 * log dflt msg for cuda err @e plus custom string @errStr
 * return POSIX - ZERO on no error occurred
 */
inline int cudaErr(cudaError_t e,const char* errStr){
	if(e){
		ERRPRINTS("%s\t%s\n",errStr,cudaGetErrorString(e));
		checkCudaErrors(e); //TODO code=700(cudaErrorIllegalAddress) "e" ??
	}
	return e; //!= 0;
}

///SPMAT UPLOAD UTILS
const cudaMemcpyKind dirUp		= cudaMemcpyHostToDevice;
const cudaMemcpyKind dirDown	= cudaMemcpyDeviceToHost;
//#define dirUp	 cudaMemcpyHostToDevice;
//#define dirDown	 cudaMemcpyDeviceToHost;

/*
 * copy CSR stored matrix in @m (Host Pntr) inside GPU pointer in @dst
 */
int spMatCpyCSR(spmat* m,spmat* dst);
///
/*
 * copy ELL stored matrix in @m (Host Pntr) inside GPU pointer in @dst
 * pitches of JA and AS ebbedded in specific fields of *dst
 */
int spMatCpyELL(spmat* m,spmat* dst);
//version with the malloc/memcpy with pitching
int spMatCpyELLNNPitched(spmat* m,spmat* dst);

inline int cudaFreeSpmat(spmat* m){
	int out = cudaFree(m->AS);
	out    |= cudaFree(m->JA);
	out    |= cudaFree(m->IRP);	///TODO 0MEMSET REQUIRED! for ELL
	#ifdef ROWLENS
	out    |= cudaFree(m->RL);
	#endif
	return out;
}

///VAR MISC

/*
 * fill @vShrd vector in shared mem with corrisponding values in @v 
 * using scheduled thread in a 1D grid of 1D blocks
 * 	-> use every thread available, reiterating if the whole grid is not enough
 */
inline __device__ void fillSharedVector_block1D_grid1D(double* v, double* vShrd, ulong len){ //TODO inline
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

//REDUCE
inline __device__ double reduceWarpRegs(double x) {
   for (int off=16; off>0; off/=2) {
      x += __shfl_down_sync(0xffffffff, x, off, 32);
   }
   return x;
}

#endif //UTILSCUDA_H
