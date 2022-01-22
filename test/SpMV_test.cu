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
#include "parser.h"
#include "utils.h"
#include "macros.h"
#include "ompChunksDivide.h"
#include "ompGetICV.h"  //ICV - RUNTIME information audit auxs
#ifdef CBLAS_TESTS
//serial dense computation with CBLAS
double* SGEMVCBLAS(spmat* mat, double* inVect);
#endif

CHUNKS_DISTR    chunksFair,chunksFairFolded,chunksNOOP;

#ifdef __CUDACC__
}						//extern "C"
//__constant__ struct config ConfD; TODO???
#include <cuda_runtime.h>  // For CUDA runtime API
//#include <helper_cuda.h>   // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers
#include "cudaUtils.h"
#endif
///inline exports
//SPMV_CHUNKS_DISTR spmvChunksFair; 
//spmat* allocSpMatrix(ulong rows, ulong cols);
//int allocSpMatrixInternal(ulong rows, ulong cols, spmat* mat);
//void freeSpmatInternal(spmat* mat);
//void freeSpmat(spmat* mat);

//global vars	->	audit
double Start,End,Elapsed,ElapsedInternal;
CHUNKS_DISTR_INTERF chunkDistrbFunc=&chunksFair; //Folded;
CONFIG Conf = {
    .gridRows = 8,
    .gridCols = 8,
};

////wrap result check and stats gather of SpMV implementation func at f
static inline int testSpMVImpl(SPMV_INTERF f,spmat* mat,double* vector,
  double* outV, double* oracleOut){
    //elapsed stats aux vars
    double times[AVG_TIMES_ITERATION],  timesInteral[AVG_TIMES_ITERATION];
    double elapsedStats[2],  elapsedInternalStats[2], start,end;
    for (uint i=0;  i<AVG_TIMES_ITERATION; i++){
        start = omp_get_wtime();
        if (f(mat,vector,&Conf,outV)){
            ERRPRINTS("compute func at:%p failed...\n",f);
            return EXIT_FAILURE;
        }
        end = omp_get_wtime();
        if (doubleVectorsDiff(oracleOut,outV,mat->M,NULL)) return EXIT_FAILURE;
        times[i]        = end - start;
        timesInteral[i] = ElapsedInternal;
        ElapsedInternal = Elapsed = 0;
    }
    statsAvgVar(times,AVG_TIMES_ITERATION,elapsedStats);
    statsAvgVar(timesInteral,AVG_TIMES_ITERATION,elapsedInternalStats);
    printf("timeAvg:%le timeVar:%le\ttimeInternalAvg:%le timeInternalVar:%le \n",
      elapsedStats[0],elapsedStats[1],elapsedInternalStats[0],elapsedInternalStats[1]);
    return EXIT_SUCCESS;
}
#ifdef __CUDACC__
static inline int testSpMVImplCuda(SPMV_CUDA_INTERF f,spmat* dMat, spmat* hmat,
  double* dVect, double* dOutV, double* hOutV, double* oracleOut){
    //elapsed stats aux vars
    double times[AVG_TIMES_ITERATION],  timesInteral[AVG_TIMES_ITERATION];
    double elapsedStats[2],  elapsedInternalStats[2];
    for (uint i=0;  i<AVG_TIMES_ITERATION; i++){
		StopWatchInterface* timer = 0; sdkCreateTimer(&timer);
		timer->start();

		f<<<Conf.gridSize,Conf.blockSize>>>(dMat,dVect,Conf,dOutV);
		DEBUG{ 
			if( cudaPeekAtLastError() )	ERRPRINT("kernelLaunchErrd\n");
			if(cudaErr(cudaPeekAtLastError(),"kernelLaunchErrd"))
				return EXIT_FAILURE;
		}
		cudaError_t dSyncOut = cudaDeviceSynchronize();
		timer->stop();
		if(cudaErr(dSyncOut,"cudaDeviceSynchronize"))
			return EXIT_FAILURE;
		Elapsed = ElapsedInternal = timer->getTime() / 1e3;

		//copy back result
		//memset(hOutV,0,sizeof(*hOutV) * hmat->M);	  //clean from before
		if(cudaErr(cudaMemcpy(hOutV,dOutV,sizeof(*hOutV)*hmat->M,dirDown),
		  "CUDA downloadResoult"))
			return EXIT_FAILURE;
        //TODO !!! if (doubleVectorsDiff(oracleOut,hOutV,hmat->M,NULL)) return EXIT_FAILURE;
        doubleVectorsDiff(oracleOut,hOutV,hmat->M,NULL);
        times[i]        = Elapsed;
        timesInteral[i] = ElapsedInternal;

		timer->reset();
        ElapsedInternal = Elapsed = 0;
    }
    statsAvgVar(times,AVG_TIMES_ITERATION,elapsedStats);
    statsAvgVar(timesInteral,AVG_TIMES_ITERATION,elapsedInternalStats);
    printf("timeAvg:%le timeVar:%le\ttimeInternalAvg:%le timeInternalVar:%le \n",
      elapsedStats[0],elapsedStats[1],elapsedInternalStats[0],elapsedInternalStats[1]);
    return EXIT_SUCCESS;
}
#endif

#define TESTTESTS   "TESTTESTS"
#define RNDVECT     "RNDVECT"
#define HELP "usage: MatrixMarket_sparse_matrix_COO[.COMPRESS_EXT]," \
    " vectorFile || "RNDVECT", ["TESTTESTS" (Requires#-DCBLAS_TESTS) ]\n"
int main(int argc, char** argv){
    int out=EXIT_FAILURE;
    if (init_urndfd())  return out;
    if (argc < 3 )  {ERRPRINT(HELP); return out;}
    
    double *vector = NULL, *outV = NULL, *oracleOut=NULL;
    ulong vectSize;
    spmat *mat = NULL, *matCSR = NULL, *matELL = NULL;
	#ifdef __CUDACC__ 
	double* dVect 		= NULL;
	double* dOutV 		= NULL;
	spmat* dMat   		= NULL;
	spmat  dMatCopy;	//host copy of CUDA Device's spmat for internal free
	spmat* matELL_t   	= NULL; //transposed matrix ptr
	spmat* matELLTrgt 	= NULL; //cudaELL,aux trgt host matPntr to use
	memset(&dMatCopy,0,sizeof(dMatCopy)); //TODO DEBUG
	ulong rowsComputeNum = 0;
	#endif

    ////parse sparse matrix and dense vector
    //extract compressed matrix
    char* trgtMatrix = TMP_EXTRACTED_MARTIX;
    if (extractInTmpFS(argv[1],TMP_EXTRACTED_MARTIX) < 0)   trgtMatrix = argv[1];
    //parse for CSR implementations
    if (!(matCSR = MMtoCSR(trgtMatrix)))	return out;
    if (!(matELL = MMtoELL(trgtMatrix)))    ERRPRINTS("ELL not feasible for %s:(",argv[1]);
		//abort later in other (reallocation if needed)
	mat = matCSR;
    VERBOSE printf("parsed matrix %lu x %lu -- %lu NNZ \n",mat->M,mat->N,mat->NZ);
    vectSize = mat->N;  //size for GEMV
    int maxThreads = omp_get_max_threads();
    ////get the vector
    if (!(strncmp(argv[2],RNDVECT,strlen(RNDVECT)))){ //generate a random vector
        if (!(vector = (typeof(vector)) malloc(vectSize * sizeof(*vector)))){
            ERRPRINT("rnd vector malloc failed\n");
            goto _free;
        }
        if (fillRndVector(vectSize,vector)){
            ERRPRINT("fillRndVector errd\n");
            goto _free;
        }
        DEBUG{ 
            if (writeDoubleVectorAsStr(RNDVECTORDUMP,vector,vectSize)) 
                ERRPRINT("couldn't write rnd vector\n");
            if (writeDoubleVector(RNDVECTORDUMPRAW,vector,vectSize)) 
                ERRPRINT("couldn't write rnd RAW vector\n");
        }
    } else{ //read vector from the given file
        if (!(vector = readDoubleVector(argv[2],&vectSize))){
            fprintf(stderr,"err during readDoubleVector at:%s\n",argv[2]);
            goto _free;
        }
        CONSISTENCY_CHECKS{
            if (vectSize != mat->N){
                ERRPRINT("vector not compatible with sparse matrix\n");
                goto _free;
            }
        }
        VERBOSE printf("parsed vector of %lu NNZ \n",vectSize);
    }
    if (!(outV = (typeof(outV)) malloc( mat->M * sizeof(*outV) ))){
        ERRPRINT("outV malloc errd\n");
        goto _free;
    }
    
    DEBUGPRINT {
        printf("sparse matrix:\n");printSparseMatrix(mat,TRUE);
        printf("vector:\n");printVector(vector,vectSize);
    }
    //// get TEST OUTPUT
#ifdef CBLAS_TESTS
    //SERIAL - LAPACK.CBLAS COMPUTATION
    if (!(oracleOut = SGEMVCBLAS(mat,vector))){
        ERRPRINT("LAPACK.CBLAS SERIAL, DENSE GEMV TEST FAILED!!\n");
        goto _free;
    }
    if ( argc == 4 && !strncmp(argv[3],TESTTESTS,strlen(TESTTESTS)) ){
        printf("testing the tests:\n matching dense CBLAS test with serial impl \n");
        if (!(outV = (typeof(outV)) malloc(mat->M * sizeof(*oracleOut)))){
            ERRPRINT("outV malloc errd for serial implementation\n");
            goto _free;
        }
        sgemvSerial(mat,vector,&Conf,outV);
        out = doubleVectorsDiff(oracleOut,outV,mat->M,NULL);    
        goto _free;
    }
#else
    //SERIAL - EASY IMPLEMENTATION
    if (!(oracleOut = (typeof(oracleOut)) malloc(mat->M * sizeof(*oracleOut)))){
        ERRPRINT("oracleOut malloc errd for serial implementation\n");
        goto _free;
    }
    sgemvSerial(mat,vector,&Conf,oracleOut);
#endif
    DEBUG{
        if (writeDoubleVector(OUTVECTORDUMPRAW,oracleOut,vectSize))
            ERRPRINT("couldn't dump out\n");
        if (writeDoubleVectorAsStr(OUTVECTORDUMP,oracleOut,vectSize))
            ERRPRINT("couldn't dump out RAW\n");
    }
    if (!getConfig(&Conf)){
        VERBOSE printf("configuration changed from env\n");
    }
    printf("SpMV_OMP_test.c\tAVG_TIMES_ITERATION:%d\t"
      "sparse matrix: %lux%lu-%luNNZ-%ld=MAX_ROW_NZ - grid: %ux%u\n",
      AVG_TIMES_ITERATION,mat->M,mat->N,mat->NZ,matELL?matELL->MAX_ROW_NZ:-1,
	  Conf.gridRows,Conf.gridCols);
    //extra configuration
    Conf.threadNum = (uint) maxThreads;
    DEBUG   printf("omp_get_max_threads:\t%d\n",maxThreads); 
    /*
     * get exported schedule configuration, 
     * if schedule != static -> dyn like -> set a chunk division function before omp for
     */
    int schedKind_chunk_monotonic[3];
    ompGetRuntimeSchedule(schedKind_chunk_monotonic);
    Conf.chunkDistrbFunc = 		(void*) chunksNOOP; 
    if (schedKind_chunk_monotonic[0] != omp_sched_static)
        Conf.chunkDistrbFunc = 	(void*) chunkDistrbFunc;
    //// PARALLEL COMPUTATIONs TO CHECK
	///CUDA IMPLEMENTATIONS
	#ifdef __CUDACC__
	DEBUG  	hprintf("\n\ntesting CUDA IMPLEMENTATIONS\n");
	///copy the spMat on the device
	if(cudaErr( cudaMalloc(&dMat,sizeof(*dMat)),"dMat"))			goto _free;
	if(cudaErr( cudaMalloc(&dVect,sizeof(*dVect)*mat->N),"dVect"))	goto _free;
	if(cudaErr( cudaMalloc(&dOutV,sizeof(*dOutV)*mat->M),"dOutV"))	goto _free;
	if((out = cudaErr(cudaMemcpy(dVect,vector,sizeof(*dVect)*mat->N,dirUp),"dVectUp")))
		goto _free;
	rowsComputeNum = mat->M;
	//1thread per row->col	1D-1D threading, easiest implementation for both CSR,CUDA
	//Conf.sharedMemSize= 		sizeof(*dOutV) * mat->M;
	Conf.blockSize		= 		dim3( BLOCKS_1D );
	Conf.gridSize 		= 		dim3( INT_DIV_CEIL(rowsComputeNum,BLOCKS_1D) );
	//if((out = cudaErr(cudaMemcpy(&ConfD,&Conf,sizeof(Conf),dirUp),"dConfUp")))			goto _free;
	SPMV_CUDA_INTERF funcCuda;

	///CSR IMPLEMENTATIONs
	if(!matCSR){
    	if (!(matCSR = MMtoCSR(trgtMatrix)))	goto _free;
	}
	if(spMatCpyCSR(matCSR,dMat))	goto _free;
	//get a copy of the internal mat for later free		
	if(cudaErr(cudaMemcpy(&dMatCopy,dMat,sizeof(*dMat),dirDown),"dMatCopy"))
		goto _free;

    for (uint f=0; f<STATIC_ARR_ELEMENTS_N(SpmvCUDA_CSRFuncs); f++){
		if(f >= SpmvCUDA_CSRFuncs_WarpPerRowIdx ){
			//2D - BLOCKING->GRID 	~	 WARPIZED VERSION
			Conf.blockSize	= 	dim3( WARPSIZE, BLOCKS_2D_WARP_R );
			Conf.gridSize	= 	dim3( INT_DIV_CEIL(mat->M,BLOCKS_2D_WARP_R) );
		}
        funcCuda = SpmvCUDA_CSRFuncs[f];
        hprintsf("@computing SpMV   with func:\%u CUDA CSR at:%p\t",f,funcCuda);
        if(testSpMVImplCuda(funcCuda,dMat,matCSR,dVect,dOutV,outV,oracleOut))
			goto _free;
    }
	assert( !(cudaFreeSpmat(&dMatCopy)));
	
	///ELL IMPLEMENTATIONS
	if (!matELL){
    	if (!(matELL = MMtoELL(trgtMatrix)))   goto _free; //TODO CUDA ONLY
	}
	//first impl. uses transposed matrix for gMEM coalescing
	if(!(matELL_t = ellTranspose(matELL)))		goto _free;
	matELLTrgt = matELL_t;
	if(spMatCpyELL(matELLTrgt,dMat))			goto _free;
	//if( mat != ellMat )	freeSpmat(ellMat);	//ELL impl. coalesced ...  Trsp free 
	//get a copy of the internal mat for later free		
	if(cudaErr(cudaMemcpy(&dMatCopy,dMat,sizeof(*dMat),dirDown),"dMatCopy"))
		goto _free;
    for(uint f=0; f<STATIC_ARR_ELEMENTS_N(SpmvCUDA_ELLFuncs); f++){
		if(SpmvCUDA_ELLFuncs_NN_TraposedImpl <= f && f < SpmvCUDA_ELLFuncs_WarpPerRowIdx){	
			//NN transposed implementations
			assert( !(cudaFreeSpmat(&dMatCopy)));
			matELLTrgt = matELL;
			if(spMatCpyELL(matELLTrgt,dMat))	goto _free;
			if(cudaErr(cudaMemcpy(&dMatCopy,dMat,sizeof(*dMat),dirDown),"dMatCopy"))
				goto _free;
		}
		else if(f >= SpmvCUDA_ELLFuncs_WarpPerRowIdx ){ //NN transposed + warp threading
			//2D - BLOCKING->GRID 	~	 WARPIZED VERSION
			Conf.blockSize	= 	dim3( WARPSIZE, BLOCKS_2D_WARP_R );
			Conf.gridSize	= 	dim3( INT_DIV_CEIL(mat->M,BLOCKS_2D_WARP_R) );
		}

        hprintsf("@computing SpMV   with func:\%u CUDA ELL at:%p\t",f,funcCuda);
        funcCuda = SpmvCUDA_ELLFuncs[f];
        if(testSpMVImplCuda(funcCuda,dMat,matELLTrgt,dVect,dOutV,outV,oracleOut))  
			goto _free;
    }
	assert( !(cudaFreeSpmat(&dMatCopy)));
	//TODO NVCC MESS UP WITH OMP COMPILATION...??? MOVED HERE TO TEST ONLY CUDA
	out = EXIT_SUCCESS;	goto _free;
    #endif //__CUDACC__

	DEBUG  	hprintf("\n\ntesting OMP IMPLEMENTATIONS\n");
    SPMV_INTERF spmvFunc;
    for (uint f=0; f<STATIC_ARR_ELEMENTS_N(SpmvCSRFuncs); f++){
        spmvFunc = SpmvCSRFuncs[f];
        hprintsf("@computing SpMV   with func:\%u OMP CSR at:%p\t",f,spmvFunc);
        if(testSpMVImpl(spmvFunc,mat,vector,outV,oracleOut))  goto _free;
    }
	//#ifndef __CUDACC__	//TODO NVCC MESSED AROUND->MOVD useful for later tests
    freeSpmat(matCSR); matCSR = NULL;
	//#endif
    //ELL IMPLEMENTATIONS
	if (!matELL){
    	if (!(matELL = MMtoELL(trgtMatrix)))   goto _free; //TODO CUDA ONLY
	}
	mat = matELL;
    for (uint f=0; f<STATIC_ARR_ELEMENTS_N(SpmvELLFuncs); f++){
        spmvFunc = SpmvELLFuncs[f];
        hprintsf("@computing SpMV   with func:\%u OMP ELL at:%p\t",f,spmvFunc);
        if(testSpMVImpl(spmvFunc,mat,vector,outV,oracleOut))  goto _free;
    }

    out = EXIT_SUCCESS;
    _free:
    if (matCSR)       freeSpmat(matCSR);
    if (matELL)       freeSpmat(matELL);
    if (vector)       free(vector);
    if (outV)    free(outV);
    if (oracleOut)    free(oracleOut);
	#ifdef __CUDACC__
	_free_cuda:
	assert( !(cudaFree(dVect)));
	assert( !(cudaFree(dOutV)));
    if (matELL_t)     freeSpmat(matELL_t);
	#endif
    return out;
}
