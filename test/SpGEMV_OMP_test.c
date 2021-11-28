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

///inline export here 
//SPMV_CHUNKS_DISTR spmvChunksFair; 
spmat* allocSpMatrix(ulong rows, ulong cols);
int allocSpMatrixInternal(ulong rows, ulong cols, spmat* mat);
void freeSpmatInternal(spmat* mat);
void freeSpmat(spmat* mat);
CHUNKS_DISTR    chunksFair,chunksFairFolded,chunksNOOP;
CHUNKS_DISTR_INTERF chunkDistrbFunc=&chunksFairFolded;
#ifdef CBLAS_TESTS
#include <cblas.h>
double* SGEMVCBLAS(spmat* mat, double* inVect){
    CBLAS_LAYOUT layout=CblasRowMajor;
    CBLAS_TRANSPOSE notrans=CblasNoTrans;
    CBLAS_INT m=mat->M, n=mat->N;
    double* denseMat = CSRToDense(mat);
    if (!denseMat){
        ERRPRINT("GEMVCheckCBLAS: aux dense matrix alloc failed\n");
        return NULL;
    }
    double* oracleOut = malloc(m * sizeof(*oracleOut));
    if (!oracleOut){
        ERRPRINT("GEMVCheckCBLAS: out for serial oracle malloc failed\n");
        goto _err;
    }
    VERBOSE printf("computimg Sparse GEMV using densification over LAPACK.CBLAS\n");
    cblas_dgemv(layout,notrans,m,n,1.0,denseMat,n,inVect,1,0.0,oracleOut,1);
    
    _free:
    free(denseMat);
    return oracleOut;

    _err:
    if (denseMat)   free(denseMat);
    if (oracleOut)  free(oracleOut);
    return NULL;
}
#endif 

CONFIG Conf = {
    .gridRows = 8,
    .gridCols = 8,
};

#define TESTTESTS   "TESTTESTS"
#define RNDVECT     "RNDVECT"
#define HELP "usage: MatrixMarket_sparse_matrix_COO[.COMPRESS_EXT], vectorFile || "RNDVECT", ["TESTTESTS"]\n"
int main(int argc, char** argv){
    int out=EXIT_FAILURE;
    if (init_urndfd())  return out;
    if (argc < 3 )  {ERRPRINT(HELP); return out;}
    
    double *vector = NULL, *outVector = NULL, *oracleOut=NULL;
    ulong vectSize;
    spmat* mat = NULL; 
    ////parse sparse matrix and dense vector
    char* trgtMatrix = TMP_EXTRACTED_MARTIX;
    if (extractInTmpFS(argv[1],TMP_EXTRACTED_MARTIX) < 0)   trgtMatrix = argv[1];
    if (!(mat = MMtoCSR(trgtMatrix))){
        ERRPRINT("err during conversion MM -> CSR\n");
        return out;
    }
    ////get the vector
    if (!(strncmp(argv[2],RNDVECT,strlen(RNDVECT)))){ //generate a random vector
        vectSize = mat->N;  //size for GEMV
        if (!(vector = malloc(vectSize * sizeof(*vector)))){
            ERRPRINT("rnd vector malloc failed\n");
            goto _free;
        }
        if (fillRndVector(vectSize,vector)){
            ERRPRINT("fillRndVector errd\n");
            goto _free;
        }
    } else{ //read vector from the given file
        if (!(vector = readVector(argv[2],&vectSize))){
            fprintf(stderr,"err during readVector at:%s\n",argv[2]);
            goto _free;
        }
        CONSISTENCY_CHECKS{
            if (vectSize != mat->N){
                ERRPRINT("vector not compatible with sparse matrix\n");
                goto _free;
            }
        }
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
        out = doubleVectorsDiff(oracleOut,outVector,mat->M);    
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
    
    printf("SpGEMV_OMP_test.c\tAVG_TIMES_ITERATION:%d\t"
      "sparse matrix: %lux%lu-%luNNZ - grid: %ux%u\n",
      AVG_TIMES_ITERATION,mat->M,mat->N,mat->NZ,Conf.gridRows,Conf.gridCols);
    //// PARALLEL COMPUTATIONs TO CHECK
    //extra configuration
    int maxThreads = omp_get_max_threads();
    Conf.threadNum = (uint) maxThreads;
    /*
     * get exported schedule configuration, 
     * if chunkSize == 1 set a chunk division function before omp for
     */
    int schedKind_chunk_monotonic[3];
    ompGetRuntimeSchedule(schedKind_chunk_monotonic);
    Conf.chunkDistrbFunc = chunksNOOP; 
    if (schedKind_chunk_monotonic[1] == 1)  Conf.chunkDistrbFunc = chunkDistrbFunc;
    if (!getConfig(&Conf)){
        VERBOSE printf("configuration changed from env");
    }
    //elapsed stats aux vars
    double times[AVG_TIMES_ITERATION],  timesInteral[AVG_TIMES_ITERATION];
    double elapsedStats[2],  elapsedInternalStats[2], start,end;
    uint f;
    SPGEMV_INTERF spgemvFunc;
    for (f=0,spgemvFunc=SpgemvFuncs[f]; spgemvFunc; spgemvFunc=SpgemvFuncs[++f]){
        hprintsf("@computing SpGEMV   with func:\%u at:%p\t",f,spgemvFunc);
        for (uint i=0;  i< AVG_TIMES_ITERATION; i++){
            start = omp_get_wtime();
            if (spgemvFunc(mat,vector,&Conf,outVector)){
                ERRPRINTS("compute func number:%u failed...\n",f);
                goto _free;
            }
            end = omp_get_wtime();
            if ((out = doubleVectorsDiff(oracleOut,outVector,mat->M))) goto _free;
            times[i]        = end - start;
            timesInteral[i] = ElapsedInternal;
            ElapsedInternal = Elapsed = 0;
        }
        statsAvgVar(times,AVG_TIMES_ITERATION,elapsedStats);
        statsAvgVar(timesInteral,AVG_TIMES_ITERATION,elapsedInternalStats);
        printf("timeAvg:%le timeVar:%le\ttimeInternalAvg:%le timeInternalVar:%le \n",
          elapsedStats[0],elapsedStats[1],elapsedInternalStats[0],elapsedInternalStats[1]);
    }

    _free:
    if (mat)          freeSpmat(mat);
    if (vector)       free(vector);
    if (outVector)    free(outVector);
    if (oracleOut)    free(oracleOut);
    return out;
}
