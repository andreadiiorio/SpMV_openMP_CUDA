#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "sparseMatrix.h"
#include "SpGEMV.h"
#include "parser.h"
#include "utils.h"
#include "macros.h"
#include "ompChunksDivide.h"
#include "ompGetICV.h"  //ICV - RUNTIME information audit auxs
#ifdef CBLAS_TESTS
//serial dense computation with CBLAS
double* SGEMVCBLAS(spmat* mat, double* inVect);
#endif

///inline exports
//SPMV_CHUNKS_DISTR spmvChunksFair; 
//spmat* allocSpMatrix(ulong rows, ulong cols);
//int allocSpMatrixInternal(ulong rows, ulong cols, spmat* mat);
//void freeSpmatInternal(spmat* mat);
//void freeSpmat(spmat* mat);

CHUNKS_DISTR    chunksFair,chunksFairFolded,chunksNOOP;

//global vars	->	audit
double Start,End,Elapsed,ElapsedInternal;
CHUNKS_DISTR_INTERF chunkDistrbFunc=&chunksFair; //Folded;
CONFIG Conf = {
    .gridRows = 8,
    .gridCols = 8,
};

//wrap result check and stats gather of SpMV implementation func at f
static inline int testSpMVImpl(SPGEMV_INTERF f,spmat* mat,double* vector,
  double* outVector, double* oracleOut){
    //elapsed stats aux vars
    double times[AVG_TIMES_ITERATION],  timesInteral[AVG_TIMES_ITERATION];
    double elapsedStats[2],  elapsedInternalStats[2], start,end;
    for (uint i=0;  i<AVG_TIMES_ITERATION; i++){
        start = omp_get_wtime();
        if (f(mat,vector,&Conf,outVector)){
            ERRPRINTS("compute func at:%p failed...\n",f);
            return EXIT_FAILURE;
        }
        end = omp_get_wtime();
        if (doubleVectorsDiff(oracleOut,outVector,mat->M,NULL)) return EXIT_FAILURE;
        times[i]        = end - start;
        timesInteral[i] = ElapsedInternal;
        ElapsedInternal = Elapsed = 0;
    }
    statsAvgVar(times,AVG_TIMES_ITERATION,elapsedStats);
    statsAvgVar(timesInteral,AVG_TIMES_ITERATION,elapsedInternalStats);
    printf("timeAvg:%le timeVar:%le\ttimeInternalAvg:%le timeInternalVar:%le \n",
      elapsedStats[0],elapsedStats[1],elapsedInternalStats[0],elapsedInternalStats[1]);
    return 0;
}

#define TESTTESTS   "TESTTESTS"
#define RNDVECT     "RNDVECT"
#define HELP "usage: MatrixMarket_sparse_matrix_COO[.COMPRESS_EXT]," \
    " vectorFile || "RNDVECT", ["TESTTESTS" (Requires#-DCBLAS_TESTS) ]\n"
int main(int argc, char** argv){
    int out=EXIT_FAILURE;
    if (init_urndfd())  return out;
    if (argc < 3 )  {ERRPRINT(HELP); return out;}
    
    double *vector = NULL, *outVector = NULL, *oracleOut=NULL;
    ulong vectSize;
    spmat* mat = NULL; 
    ////parse sparse matrix and dense vector
    //extract compressed matrix
    char* trgtMatrix = TMP_EXTRACTED_MARTIX;
    if (extractInTmpFS(argv[1],TMP_EXTRACTED_MARTIX) < 0)   trgtMatrix = argv[1];
    //parse for CSR implementations
    if (!(mat = MMtoCSR(trgtMatrix))){
        ERRPRINT("err during conversion MM -> CSR\n");
        return out;
    }
    VERBOSE printf("parsed matrix %lu x %lu -- %lu NNZ \n",mat->M,mat->N,mat->NZ);
    vectSize = mat->N;  //size for GEMV
    ////get the vector
    if (!(strncmp(argv[2],RNDVECT,strlen(RNDVECT)))){ //generate a random vector
        if (!(vector = malloc(vectSize * sizeof(*vector)))){
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
    if (!(outVector = malloc( mat->M * sizeof(*outVector) ))){
        ERRPRINT("outVector malloc errd\n");
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
        if (!(outVector = malloc(mat->M * sizeof(*oracleOut)))){
            ERRPRINT("outVector malloc errd for serial implementation\n");
            goto _free;
        }
        sgemvSerial(mat,vector,&Conf,outVector);
        out = doubleVectorsDiff(oracleOut,outVector,mat->M,NULL);    
        goto _free;
    }
#else
    //SERIAL - EASY IMPLEMENTATION
    if (!(oracleOut = malloc(mat->M * sizeof(*oracleOut)))){
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
    printf("SpGEMV_OMP_test.c\tAVG_TIMES_ITERATION:%d\t"
      "sparse matrix: %lux%lu-%luNNZ - grid: %ux%u\n",
      AVG_TIMES_ITERATION,mat->M,mat->N,mat->NZ,Conf.gridRows,Conf.gridCols);
    //extra configuration
    int maxThreads = omp_get_max_threads();
    Conf.threadNum = (uint) maxThreads;
    DEBUG   printf("omp_get_max_threads:\t%d\n",maxThreads); 
    /*
     * get exported schedule configuration, 
     * if schedule != static -> dyn like -> set a chunk division function before omp for
     */
    int schedKind_chunk_monotonic[3];
    ompGetRuntimeSchedule(schedKind_chunk_monotonic);
    Conf.chunkDistrbFunc = chunksNOOP; 
    if (schedKind_chunk_monotonic[0] != omp_sched_static)
        Conf.chunkDistrbFunc = chunkDistrbFunc;
    //// PARALLEL COMPUTATIONs TO CHECK
    SPGEMV_INTERF spgemvFunc;
    for (uint f=0; f<STATIC_ARR_ELEMENTS_N(SpgemvCSRFuncs); f++){
        spgemvFunc = SpgemvCSRFuncs[f];
        hprintsf("@computing SpGEMV   with func:\%u CSR at:%p\t",f,spgemvFunc);
        if(testSpMVImpl(spgemvFunc,mat,vector,outVector,oracleOut))  goto _free;
    }
    //ELL IMPLEMENTATIONS
    freeSpmat(mat);
    if (!(mat = MMtoELL(trgtMatrix)))   goto _free;
    for (uint f=0; f<STATIC_ARR_ELEMENTS_N(SpgemvELLFuncs); f++){
        spgemvFunc = SpgemvELLFuncs[f];
        hprintsf("@computing SpGEMV   with func:\%u ELL at:%p\t",f,spgemvFunc);
        if(testSpMVImpl(spgemvFunc,mat,vector,outVector,oracleOut))  goto _free;
    }
    
    out = EXIT_SUCCESS;
    _free:
    if (mat)          freeSpmat(mat);
    if (vector)       free(vector);
    if (outVector)    free(outVector);
    if (oracleOut)    free(oracleOut);
    return out;
}
