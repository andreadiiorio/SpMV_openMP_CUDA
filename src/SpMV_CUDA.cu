/*
 * 	dev:	Andrea Di Iorio
 * 	CUDA implementations of SpMV with CSR and ELL spMatrix rappresetantations
 */

extern "C" {
#include "sparseMatrix.h"
#include "SpMV.h"
}
#include "cudaUtils.h"

///inline exports
//__device__ double reduceWarpRegs(double x);


//1 threadPerRow
__global__ void cudaSpMVRowsCSR(spmat* m,double* v,CONFIG cfg,double* outV){
	uint tRow = threadIdx.x + blockIdx.x*blockDim.x;	//1D block,1D grid
	if (tRow < m->M){
		ulong eIdx,rLen;
		#ifdef ROWLENS
		rLen = m->RL[tRow];
		#else
		rLen = m->IRP[tRow+1] - m->IRP[tRow];
		#endif
		double outEntry = 0;	//thread's output entry in outV
		for(ulong c=0; c<rLen; c++){
			eIdx = m->IRP[tRow] + c;
			outEntry += m->AS[eIdx] * v[m->JA[eIdx]];
		}
		outV[tRow] = outEntry;
	}
}

//1 warpPerRow
__global__ void cudaSpMVWarpPerRowCSR(spmat* m,double* v,CONFIG cfg,double* outV){
	//uint tId  = threadIdx.x + blockIdx.x*blockDim.x;
	//uint tWarpIdx = tId % warpSize,tRow = tId / warpSize
	//1D GRID of 2D BLOCKS: [warpIdx,rowIdx]
	uint tWarpIdx = threadIdx.x;
	uint tRow = threadIdx.y + blockIdx.y*blockDim.y;
	if (tRow < m->M){
		ulong eIdx,rLen;
		#ifdef ROWLENS
		rLen = m->RL[tRow];
		#else
		rLen = m->IRP[tRow+1] - m->IRP[tRow];
		#endif
		double outEntry = 0;	//thread's output entry in outV
		for(ulong c=tWarpIdx; c<rLen; c+=warpSize){
			eIdx = m->IRP[tRow] + c;
			outEntry += m->AS[eIdx] * v[m->JA[eIdx]];
		}
		outEntry = reduceWarpRegs(outEntry);	//in-reg&inter warp tree reduction
		if (!tWarpIdx)	outV[tRow] = outEntry;	//0th thread will see the reducted out
	}
}

///ELL METHODS

//1thread per m.row -> m^T.col (coalescingUp)
//EXPECTED MATRIX @m TRANSPOSED for gMem coealscing
__global__ void cudaSpMVRowsELL(spmat* m,double* v,CONFIG cfg,double* outV){
	uint tId = threadIdx.x + blockIdx.x*blockDim.x;	//1D block,1D grid
	if (tId < m->MAX_ROW_NZ){
		double outEntry = 0;
		ulong cLen = m->M;
		ulong pitchR = m->pitchAS;	//TODO same for JA if same element sizes!
		//for each:	i=col => transpose => row
		for(ulong i=0,asIdx=IDX2D(i,tId,m->pitchAS),jaIdx=IDX2D(i,tId,m->pitchJA);
			i<cLen;	i++,asIdx+=pitchR,jaIdx+=pitchR){
			outEntry += m->AS[asIdx] * v[m->JA[jaIdx]];
		}		
		outV[tId] = outEntry;
	}
}

//1thread per m.row -> m^T.col (coalescingUp)
__global__ void cudaSpMVRowsELLNNTransposed(spmat* m,double* v,CONFIG cfg,double* outV){
	uint tRow = threadIdx.x + blockIdx.x*blockDim.x;	//1D block,1D grid
	if (tRow < m->M){
		double outEntry = 0;
		ulong rLen;
		#ifdef ROWLENS
		rLen = m->RL[tRow];
		#else
		rLen = m->MAX_ROW_NZ;
		#endif
		for(ulong c=0,asIdx=IDX2D(tRow,c,m->pitchAS),jaIdx=IDX2D(tRow,c,m->pitchJA);
			c<rLen;	 c++,asIdx++,jaIdx++){
				outEntry += m->AS[asIdx] * v[m->JA[jaIdx]];
		}		
		outV[tRow] = outEntry;
	}
}
__global__ void cudaSpMVWarpsPerRowELLNTrasposed(spmat* m,double* v,CONFIG cfg,double* outV){
	//1D GRID of 2D BLOCKS: [warpIdx,rowIdx]
	uint tWarpIdx = threadIdx.x;
	uint tRow = threadIdx.y + blockIdx.y*blockDim.y;
	if (tRow < m->M){
		double outEntry = 0;
		ulong rLen;
		#ifdef ROWLENS
		rLen = m->RL[tRow];
		#else
		rLen = m->MAX_ROW_NZ;
		#endif
		for(ulong c=tWarpIdx,asIdx=IDX2D(tRow,c,m->pitchAS),jaIdx=IDX2D(tRow,c,m->pitchJA);
			c<rLen;	 c+=warpSize,asIdx+=warpSize,jaIdx+=warpSize){
				outEntry += m->AS[asIdx] * v[m->JA[jaIdx]];
		}		
		outEntry = reduceWarpRegs(outEntry);	//in-reg&inter warp tree reduction
		if (!tWarpIdx)	outV[tRow] = outEntry;	//0th thread will see the reducted out
	}
}

template <uint warpsNum> __global__ void cudaSpMVWarpsPerRowELLNTrasposed(spmat* m,double* v,CONFIG cfg,double* outV){
	//1D GRID of 2D BLOCKS: [groupWarpsIdx,rowIdx]
	uint tWarpsIdx = threadIdx.x;
	uint tRow = threadIdx.y + blockIdx.y*blockDim.y;
	if (tRow< m->M){
		extern __shared__ double outEntries[]; //per thread in warpsGroup
		ulong rLen;
		#ifdef ROWLENS
		rLen = m->RL[tRow];
		#else
		rLen = m->MAX_ROW_NZ;
		#endif
		//each thread in the warpsGroup multiply an elements
		for(ulong c=tWarpsIdx, asIdx=IDX2D(tRow,c,m->pitchAS),jaIdx=IDX2D(tRow,c,m->pitchJA);
			c<rLen;		c+=blockDim.x,asIdx+=blockDim.x,jaIdx+=blockDim.x){
				outEntries[tWarpsIdx] += m->AS[asIdx] * v[m->JA[jaIdx]];
		}
		//TODO outEntries = gpu_reduce(
		//TODO if (!tWarpsIdx)	outV[tRow] = outEntry;	//0th thread will see the reducted out
	}
}
