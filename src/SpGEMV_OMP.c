#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "SpGEMV.h"
#include "parser.h"
#include "utils.h"
#include "macros.h"
#include "sparseMatrix.h"

//TODO MOVE TO HEDER
//COMPUTE MODES STRINGS

#define _ROWS            "ROWS"
#define _SORTED_ROWS     "SORTED_ROWS"
#define _TILES           "TILES"

typedef enum {
    ROWS,
    SORTED_ROWS,
    TILES
} COMPUTE_MODE;

typedef struct{
    ushort gridRows;
    ushort gridCols;
    //TODO FULL CONFIG DOCCED HERE
    int threadNum;  //thread num to use in an OMP parallel region ...
} CONFIG;
CONFIG Conf = {
    .gridRows = 8,
};

//compute function interface and its pointer definitions
typedef int (COMPUTEFUNC) (spmat*,double*,CONFIG*,double*);
typedef int (*COMPUTEFUNC_INTERF) (spmat*,double*,CONFIG*,double*);
COMPUTEFUNC spgemvRowsBasic;

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
    COMPUTEFUNC_INTERF computeFunc;
    switch (cmode){
        case ROWS:         computeFunc=&spgemvRowsBasic;break;
        case SORTED_ROWS:  printf("s");break;
        case TILES:  printf("t");break;
    }
    
    double *vector = NULL, *outVector = NULL;
    uint vectSize;
    spmat* mat = NULL; 
    ////parse sparse matrix and dense vector
    if (!(mat = MMtoCSR(argv[1]))){
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
    
    DEBUG {
        printf("sparse matrix:\n");printSparseMatrix(mat,TRUE);
        printf("vector:\n");printVector(vector,vectSize);
    }
  
    //// PARALLEL COMPUTATION
    int maxThreads = omp_get_max_threads();
    if ((out = computeFunc(mat,vector,&Conf,outVector))){
        ERRPRINT("compute function selected failed...\n"); goto _free;
    }
    printf("outVector:\n");printVector(outVector,vectSize);


    _free:
    if (mat)          free(mat);
    if (vector)       free(vector);
    if (outVector)    free(outVector);
    return out;
}

//TODO CATALOG COMPUTE FUNC IN SEP. FILES IF NEEDED
/*
 * basic spgemv with row partitioning, 1 row per thread in consecutive order
 * statically assigned to threads
 */
int spgemvRowsBasic(spmat* mat, double* vect, CONFIG* cfg, double* outVect){
    int out = EXIT_FAILURE;
    
    double acc,end,start;
    start = omp_get_wtime();

    uint i,j,k;
    omp_set_num_threads(cfg->threadNum);
    #pragma omp parallel for schedule(static) shared(vect,mat) private(i, j, k)
    for (i=0;i<mat->M;i++){
        acc = 0;
        for (j=mat->IRP[i]; j<mat->IRP[i+1]; j++){
           k = mat->JA[j];
           acc += mat->AS[j] * vect[k];
        }
        outVect[i] = acc;
    }
 
    end = omp_get_wtime();
    VERBOSE printf("spgemvRowsBasic with %u x %u CSR sp.Mat, elapsed %lf\n",
        mat->M,mat->N,end-start);
    out = EXIT_SUCCESS;
    return out;
}
