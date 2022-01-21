/*	DEV:	ANDREA DI IORIO
 *  SpGEMV  multi implementation simple cli interface
 *  wrapped CUDA specific code in nvcc exported macro __CUDACC__
 * 		-> easily excludable from CPU-ONLY compilation
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#ifdef __CUDACC__
	#pragma message("compiling with NVCC" )
	extern "C" {		//remainig function declared with C linkage
#else
	#include "cuda_c_gcc_mock.h"
#endif
#include "sparseMatrix.h"
#include "SpMV.h"
#include "ompChunksDivide.h"
#include "parser.h"
#include "utils.h"
#include "macros.h"

CHUNKS_DISTR    	chunksFair,chunksFairFolded,chunksNOOP;
CHUNKS_DISTR_INTERF chunkDistrbFunc=&chunksFairFolded;

#ifdef __CUDACC__
}
//__constant__ struct config ConfD; TODO???
#include <cuda_runtime.h>  // For CUDA runtime API
//#include <helper_cuda.h>   // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers
#include "cudaUtils.h"
#endif

//inline funcs
//void freeSpmatInternal(spmat* mat);

//global vars	->	audit
double Start,End,Elapsed,ElapsedInternal;
CONFIG Conf = {
    .gridRows = 8,
    .gridCols = 8,
};

#define RNDVECT "RNDVECT"
#define COMPUTE_MODES_OMP  "\n\tOMP:\t"CSR_ROWS","CSR_ROWS_GROUPS","CSR_TILES","CSR_TILES_ALLOCD\
	","ELL_ROWS","ELL_ROWS_GROUPS","ELL_TILES
#define COMPUTE_MODES_CUDA "\n\tCUDA:\t"CUDA_CSR_ROWS","CUDA_CSR_ROWS_WARP","CUDA_ELL_ROWS","CUDA_ELL_ROWS_WARP_NT
#define HELP "usage: MatrixMarket_sparse_matrix_COO, vectorFile || "RNDVECT", [COMPUTE MODE:"\
	COMPUTE_MODES_OMP COMPUTE_MODES_CUDA

int main(int argc, char** argv){
    int out=EXIT_FAILURE;
    if (init_urndfd())  return out;
    if (argc < 3 )  {ERRPRINT(HELP); return out;}
    ///init
    double *vector = NULL, *outV = NULL, start,end,elapsed;
    ulong vectSize;
    spmat* mat = NULL; 
	SPMV_INTERF 	   func			 = NULL;
	int toCSR = 1;	//select a CSR method
	#ifdef __CUDACC__ 
    SPMV_CUDA_INTERF funcCuda		 = NULL;
	double* dVect = NULL;
	double* dOutV = NULL;
	spmat* dMat   = NULL;
	spmat dMatCopy;	//host copy of CUDA Device's spmat for internal free
	memset(&dMatCopy,0,sizeof(dMatCopy)); //TODO DEBUG
	#endif
    //extra configuration
    int maxThreads = omp_get_max_threads();
    Conf.threadNum = (uint) maxThreads;
    /*
     * get exported schedule configuration, 
     * if chunkSize == 1 set a chunk division function before omp for
     */
    int schedKind_chunk_monotonic[3];
    ompGetRuntimeSchedule(schedKind_chunk_monotonic);
    Conf.chunkDistrbFunc = (typeof(Conf.chunkDistrbFunc)) chunksNOOP; 
    if (schedKind_chunk_monotonic[1] == 1)  
		Conf.chunkDistrbFunc = (typeof(Conf.chunkDistrbFunc)) chunkDistrbFunc;
    if (!getConfig(&Conf)){
        VERBOSE printf("configuration changed from env");
    }
    ///set compute mode
    COMPUTE_MODE cmode = _CSR_ROWS;
    if (argc > 3 ){ //parse from argv
		char* cm = argv[3];
        if ( strEqual(cm,CSR_ROWS) )				 cmode = _CSR_ROWS;
        else if ( strEqual(cm,CSR_ROWS_GROUPS) )	 cmode = _CSR_ROWS_GROUPS;
        else if ( strEqual(cm,CSR_TILES) )			 cmode = _CSR_TILES;
        else if ( strEqual(cm,CSR_TILES_ALLOCD) )	 cmode = _CSR_TILES_ALLOCD;
        else if ( strEqual(cm,ELL_ROWS) )			 cmode = _ELL_ROWS;
        else if ( strEqual(cm,ELL_ROWS_GROUPS) )	 cmode = _ELL_ROWS_GROUPS;
        else if ( strEqual(cm,ELL_TILES) )			 cmode = _ELL_TILES;

        else if ( strEqual(cm,CUDA_CSR_ROWS_WARP) )	 cmode = _CUDA_CSR_ROWS_WARP;
        else if ( strEqual(cm,CUDA_CSR_ROWS) )		 cmode = _CUDA_CSR_ROWS;
        else if ( strEqual(cm,CUDA_ELL_ROWS_WARP_NT))cmode = _CUDA_ELL_ROWS_WARP_NT;
        else if ( strEqual(cm,CUDA_ELL_ROWS_WARP) )	 cmode = _CUDA_ELL_ROWS_WARP;
        else if ( strEqual(cm,CUDA_ELL_ROWS) )		 cmode = _CUDA_ELL_ROWS;
		
		else{ ERRPRINT("INVALID COMPUTE_MODE ARGV[3] GIVEN\n" HELP);exit(1); }
    }
   
    switch (cmode){
		case _CSR_ROWS_GROUPS:	func = &spmvRowsBlocksCSR;			break;
		case _CSR_TILES:		func = &spmvTilesCSR;				break;
		case _CSR_TILES_ALLOCD:	func = &spmvTilesAllocdCSR;			break;
		case _CSR_ROWS:			func = &spmvRowsBasicCSR;			break;
		case _ELL_ROWS:			func = &spmvRowsBasicELL; toCSR=0;	break;
		case _ELL_ROWS_GROUPS:	func = &spmvRowsBlocksELL;toCSR=0;	break;
		case _ELL_TILES:		func = &spmvTilesELL;		toCSR=0;	break;
		#ifdef __CUDACC__		///CUDA IMPLEMENTATIONS
		case _CUDA_CSR_ROWS:		funcCuda = &cudaSpMVRowsCSR;		break;
		case _CUDA_ELL_ROWS:		funcCuda = &cudaSpMVRowsELL;toCSR=0;break;
		case _CUDA_CSR_ROWS_WARP:	funcCuda = &cudaSpMVWarpPerRowCSR;	break;
		case _CUDA_ELL_ROWS_WARP_NT:funcCuda = &cudaSpMVWarpsPerRowELLNTrasposed;toCSR=0;break;
		//case _CUDA_ELL_ROWS_WARP:	funcCuda = &cudaSpMVWarpPerRowELL;toCSR=0; 	break;
		#endif
		//TODO ERR ON NOT SUPPORTED OMP MODES ? 
    }
    ////parse sparse matrix and dense vector
    char* trgtMatPath = TMP_EXTRACTED_MARTIX;
    if (extractInTmpFS(argv[1],TMP_EXTRACTED_MARTIX) < 0)   trgtMatPath = argv[1];
    if ( (toCSR && !(mat = MMtoCSR(trgtMatPath))) ){
        ERRPRINT("err during parsing MatrixMarket -> CSR\n");
        return out;
    } else if ( !toCSR ){ //ELL
		if( !(mat = MMtoELL(trgtMatPath)) ){
        	ERRPRINT("err during parsing MatrixMarket -> ELL\n");
        	return out;
		}
	}
    ////get the vector
    vectSize = mat->N;
    if (!(strncmp(argv[2],RNDVECT,strlen(RNDVECT)))){ //generate a random vector
        if (!(vector = (typeof(vector)) malloc(vectSize * sizeof(*vector)))){
            ERRPRINT("rnd vector malloc failed\n");
            goto _free;
        }
        if (fillRndVector(vectSize,vector)){
            ERRPRINT("fillRndVector errd\n");
            goto _free;
        }
		VERBOSE	hprintf("dumping RNDVECT at:\t" RNDVECTORDUMP "\n");
		//if(writeDoubleVectorAsStr(RNDVECTORDUMP,outV,mat->M))	
		if(writeDoubleVector(RNDVECTORDUMP,vector,vectSize))
			ERRPRINT("RNDVECT dump err\n");
    } else{ //read vector from the given file
        if (!(vector = readDoubleVector(argv[2],&vectSize))){
        //if (!(vector = readDoubleVectorStr(argv[2],&vectSize))){
            fprintf(stderr,"err during readDoubleVectorStr at:%s\n",argv[2]);
            goto _free;
        }
        CONSISTENCY_CHECKS{
            if (vectSize != mat->N){
                ERRPRINT("vector not compatible with sparse matrix\n");
                goto _free;
            }
        }
    }
    
    if (!(outV = (typeof(outV)) malloc( mat->M * sizeof(*outV)))){
        ERRPRINT("outV malloc errd\n");
        goto _free;
    }
    
    DEBUGPRINT {
        printf("sparse matrix:\n");printSparseMatrix(mat,TRUE);
        printf("vector:\n");printVector(vector,vectSize);
    }
  
    //// PARALLEL COMPUTATION
	#ifdef __CUDACC__ 
	if( funcCuda ){		///CUDA IMPLEMENTATION
		///copy the spMat on the device
		if(cudaErr( cudaMalloc(&dMat,sizeof(*dMat)),"dMat"))			goto _free;
		if(cudaErr( cudaMalloc(&dVect,sizeof(*dVect)*mat->N),"dVect"))	goto _free;
		if(cudaErr( cudaMalloc(&dOutV,sizeof(*dOutV)*mat->M),"dOutV"))	goto _free;
		if((out = cudaErr(cudaMemcpy(dVect,vector,sizeof(*dVect)*mat->N,dirUp),"dVectUp")))
			goto _free;
		ulong rowsComputeNum = mat->M;
		//CSR IMPLEMENTATIONs
		if(cmode <= _CUDA_CSR_ROWS_WARP){
			if(spMatCpyCSR(mat,dMat))									goto _free;
		}
		//ELL IMPLEMENATIONS
		else{
			spmat* ellMat = mat;
			if(cmode < _CUDA_ELL_ROWS_WARP_NT){	//ELL impl. coalesced by transposing
				if (!(ellMat = ellTranspose(mat)))						goto _free;
			}
			if(spMatCpyELL(ellMat,dMat))								goto _free;
			if( mat != ellMat )	freeSpmat(ellMat);	//ELL impl. coalesced ...  Trsp free 
		}
		//get a copy of the internal mat for later free		
		if(cudaErr(cudaMemcpy(&dMatCopy,dMat,sizeof(*dMat),dirDown),"dMatCopy"))
			goto _free;

		//Conf.sharedMemSize= 		sizeof(*dOutV) * mat->M;
		//2D - BLOCKING->GRID 	~	 WARPIZED VERSION
		//switch(cmode){	case _CUDA_ELL_ROWS 
		Conf.blockSize	= 	dim3( WARPSIZE, BLOCKS_2D_WARP_R );
		Conf.gridSize	= 	dim3( INT_DIV_CEIL(mat->M,BLOCKS_2D_WARP_R) );
		if (cmode == _CUDA_CSR_ROWS || cmode == _CUDA_ELL_ROWS){ // NN warp methods -> 1D blocks/grid
			Conf.blockSize		= 		dim3( BLOCKS_1D );
			Conf.gridSize 		= 		dim3( INT_DIV_CEIL(rowsComputeNum,BLOCKS_1D) );
		}
		//if((out = cudaErr(cudaMemcpy(&ConfD,&Conf,sizeof(Conf),dirUp),"dConfUp")))			goto _free;

		StopWatchInterface* timer = 0; sdkCreateTimer(&timer);
		timer->start();
		//DEBUG hprintsf("inVectSize-shMem=%uKB\n",mat->M*sizeof(*vector) >> 10);
		//funcCuda<<<Conf.gridSize,Conf.blockSize,48<<10-1>>>(dMat,dVect,Conf,dOutV);
		funcCuda<<<Conf.gridSize,Conf.blockSize>>>(dMat,dVect,Conf,dOutV);
		DEBUG{ 
			if( cudaPeekAtLastError() )	ERRPRINT("kernelLaunchErrd\n");
			if((out = cudaErr(cudaPeekAtLastError(),"kernelLaunchErrd")))	goto _free;
		}
		cudaError_t dSyncOut = cudaDeviceSynchronize();
		timer->stop();
		if((out = cudaErr(dSyncOut,"cudaDeviceSynchronize")))	goto _free;
		Elapsed = ElapsedInternal = elapsed = timer->getTime();
		timer->reset();

		//copy back result
		if((out = cudaErr(cudaMemcpy(outV,dOutV,sizeof(*outV)*mat->M,dirDown),
		  "CUDA downloadResoult")))
			goto _free;
	}
	else{	
	#else	///CPU_ONLY COMPILATION
	{
	#endif
		///OMP IMPLEMENTATION SELECTED
    	start = omp_get_wtime(); 
    	if ((out = func(mat,vector,&Conf,outV))){
    	    ERRPRINT("compute function selected failed...\n"); goto _free;
    	}
		
    	end = omp_get_wtime(); elapsed = end-start;
	}
    DEBUGPRINT{
        printf("outV:\n");
        printVector(outV,vectSize);
    }
	
	VERBOSE	hprintf("dumping outV at:\t" OUTVECTORDUMP "\n");
	if(writeDoubleVector(OUTVECTORDUMPRAW,outV,mat->M) || 
	   writeDoubleVectorAsStr(OUTVECTORDUMP,outV,mat->M))	
			ERRPRINT("outV dump err\n");
    printf("cmode:%d\telapsed:\t %le elapsedInternal %le\n",
	  cmode,elapsed,ElapsedInternal);
    _free:
    if (mat) 	freeSpmat(mat);
    free(vector);
    free(outV);
	#ifdef __CUDACC__
	_free_cuda:
	assert( !(cudaFree(dVect)));
	assert( !(cudaFree(dOutV)));
	assert( !(cudaFreeSpmat(&dMatCopy)));
	#endif
    return out;
}
