#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "SpGEMV.h"
#include "ompChunksDivide.h"
#include "parser.h"
#include "utils.h"
#include "macros.h"
#include "sparseMatrix.h"

//inline funcs
CHUNKS_DISTR    chunksFair,chunksFairFolded,chunksNOOP;
CHUNKS_DISTR_INTERF chunkDistrbFunc=&chunksFairFolded;

CONFIG Conf = {
    .gridRows = 8,
    .gridCols = 8,
};

#define RNDVECT "RNDVECT"
#define HELP "usage: MatrixMarket_sparse_matrix_COO, vectorFile || "RNDVECT \
        ", [COMPUTE/PARTITION_MODE: "_ROWS","_SORTED_ROWS","_TILES" ("_ROWS")]\n"

int main(int argc, char** argv){
    int out=EXIT_FAILURE;
    if (init_urndfd())  return out;
    if (argc < 3 )  {ERRPRINT(HELP); return out;}
    ///set compute mode
    COMPUTE_MODE cmode = ROWS;
    if (argc > 3 ){ //parse from argv
        if (!(strncmp(argv[3],_SORTED_ROWS,strlen(_SORTED_ROWS))))  cmode=SORTED_ROWS;
        else if (!(strncmp(argv[3],_ROWS,strlen(_ROWS))))           cmode=ROWS;
        else if (!(strncmp(argv[3],_TILES,strlen(_TILES))))         cmode=TILES;
        else{   ERRPRINT("INVALID MODE." HELP); return out; }
    }
    SPGEMV_INTERF computeFunc;
    switch (cmode){
        case ROWS:         computeFunc=&spgemvRowsBasic;break;
        case SORTED_ROWS:  printf("s");break;
        case TILES:  printf("t");break;
    }
    
    double *vector = NULL, *outVector = NULL, start,end;
    ulong vectSize;
    spmat* mat = NULL; 
    //extra configuration
    int maxThreads = omp_get_max_threads();
    Conf.threadNum = (uint) maxThreads;
    /*
     * get exported schedule configuration, 
     * if chunkSize == 1 set a chunk division function before omp for
     */
    int schedKind_monotonic_chunk[3];
    ompGetRuntimeSchedule(schedKind_monotonic_chunk);
    Conf.chunkDistrbFunc = chunksNOOP; 
    if (schedKind_monotonic_chunk[2] == 1)  Conf.chunkDistrbFunc = chunkDistrbFunc;
    if (!getConfig(&Conf)){
        VERBOSE printf("configuration changed from env");
    }
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
    
    if (!(outVector = malloc( mat->M * sizeof(*outVector)))){
        ERRPRINT("outVector malloc errd\n");
        goto _free;
    }
    
    DEBUGPRINT {
        printf("sparse matrix:\n");printSparseMatrix(mat,TRUE);
        printf("vector:\n");printVector(vector,vectSize);
    }
  
    //// PARALLEL COMPUTATION
    start = omp_get_wtime(); 
    if ((out = computeFunc(mat,vector,&Conf,outVector))){
        ERRPRINT("compute function selected failed...\n"); goto _free;
    }
    end = omp_get_wtime(); 
    DEBUGPRINT{
        printf("outVector:\n");
        printVector(outVector,vectSize);
    }
    printf("elapsed:\t %le elapsedInternal %le\n",end-start,ElapsedInternal);

    _free:
    if (mat)          free(mat);
    if (vector)       free(vector);
    if (outVector)    free(outVector);
    return out;
}
