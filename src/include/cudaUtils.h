#ifndef UTILSCUDA_H
#define UTILSCUDA_H

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
 */
int spMatCpyELL(spmat* m,spmat* dst,size_t* pitchJA,size_t* pitchAS);

inline int cudaFreeSpmat(spmat* m){
	int out = cudaFree(m->AS);
	out    |= cudaFree(m->JA);
	out    |= cudaFree(m->IRP);	///TODO 0MEMSET REQUIRED! for ELL
	#ifdef ROWLENS
	out    |= cudaFree(m->RL);
	#endif
	return out;
}

#endif //UTILSCUDA_H
