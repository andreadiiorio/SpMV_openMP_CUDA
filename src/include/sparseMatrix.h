//sparse matrix def & aux
//TODO adapt to work on both CUDA->ELL and std CSR
#ifndef SPARSEMATRIX
#define SPARSEMATRIX 

#include "macros.h"

typedef struct{
    uint NZ,M,N;
    uint* JA;
    //CSR SPECIFIC
#ifdef ROWLENS
    uint* RL;   //row lengths
#endif
    uint* IRP;
    //CUDA SPECIFIC
    //uint MAXNZ;

    double *AS; 
} spmat; //describe a sparse matrix

////Sparse vector accumulator -- corresponding to a matrix portion
typedef struct{
    //uint    r;     //row index in the corresponding matrix
    //uint    c;     //col index in the corresponding matrix
    uint    len;   //rowLen
    double* AS;    //row nnz    values
    uint*   JA;    //row nnz    colIndexes
} SPACC; 

/*
 * return !0 if col @j idx is in row @i of sparse mat @smat
 */
inline int IS_NNZ(spmat* smat,uint i,uint j){
    int out = 0;
    for (uint x=smat->IRP[i]; x<smat->IRP[i+1] && !out; x++){
        out = (j == smat->JA[x]); 
    } 
    return out;
}
////aux functions
//free sparse matrix
inline void freeSpmatInternal(spmat* mat){
    if(mat->AS)    free(mat->AS);  
    if(mat->JA)    free(mat->JA);  
    if(mat->IRP)   free(mat->IRP);  
#ifdef ROWLENS
    if(mat->RL)    free(mat->RL);
#endif 
}
inline void freeSpmat(spmat* mat){
    freeSpmatInternal(mat);
    free(mat);
}

//free max aux structs not NULL pointed
inline void freeSpAcc(SPACC* r){ 
    if(r->AS)   free(r->AS);
    if(r->JA)   free(r->JA);
}
////alloc&init functions
//alloc&init internal structures only dependent of dimensions @rows,@cols
inline int allocSpMatrixInternal(uint rows, uint cols, spmat* mat){
    mat -> M = rows;
    mat -> N = cols;
    if (!(mat->IRP=calloc(mat->M+1,sizeof(*(mat->IRP))))){ //calloc only for 0th
        ERRPRINT("IRP calloc err\n");
        freeSpmatInternal(mat);
        return EXIT_FAILURE;
    }
#ifdef ROWLENS
    if (!(mat->RL = malloc(mat->M*sizeof(*(mat->RL))))){
        ERRPRINT("RL calloc err\n");
        freeSpmatInternal(mat);
        return EXIT_FAILURE;
    }
#endif
    return EXIT_SUCCESS;
}

//alloc a sparse matrix of @rows rows and @cols cols 
inline spmat* allocSpMatrix(uint rows, uint cols){

    spmat* mat;
    if (!(mat = calloc(1,sizeof(*mat)))) { 
        ERRPRINT("mat  calloc failed\n");
        return NULL;
    }
    if (allocSpMatrixInternal(rows,cols,mat)){
        free(mat);
        return NULL;
    }
    return mat;
}

///SPARSE MATRIX PARTITIONING
/*
 * partition CSR sparse matrix @A in @gridCols columns partitions 
 * returning an offsets matrix out[i][j] = start of jth colPartition of row i
 * subdivide @A columns in uniform cols ranges in the output 
 */
uint* colsOffsetsPartitioningUnifRanges(spmat* A,uint gridCols);

/*
 * partition CSR sparse matrix @A in @gridCols columns partitions as 
 * indipended and allocated sparse matrixes and return them
 * subdivide @A columns in uniform cols ranges in the output 
 */
spmat* colsPartitioningUnifRanges(spmat* A,uint gridCols);

/*  
    check if sparse matrixes A<->B differ up to 
    DOUBLE_DIFF_THREASH per element
*/
int spmatDiff(spmat* A, spmat* B);
////dyn alloc of spGEMM output matrix
/*
///size prediction of AB = @A * @B
inline uint SpGEMMPreAlloc(spmat* A,spmat* B){
    //TODO BETTER PREALLOC HEURISTICS HERE 
    return MAX(A->NZ,B->NZ);
}
//init a sparse matrix AB=@A * @B with a initial allocated space by an euristic
inline spmat* initSpMatrixSpGEMM(spmat* A, spmat* B){
    spmat* out;
    if (!(out = allocSpMatrix(A->M,B->N)))  return NULL;
    out -> NZ = SpGEMMPreAlloc(A,B);
    if (!(out->AS = malloc(out->NZ*sizeof(*(out->AS))))){
        ERRPRINT("initSpMatrix: out->AS malloc errd\n");
        free(out);
        return NULL;
    }
    if (!(out->JA = malloc(out->NZ*sizeof(*(out->JA))))){
        ERRPRINT("initSpMatrix: out->JA malloc errd\n");
        freeSpmat(out);
        return NULL;
    }
    return out;
}

#define REALLOC_FACTOR  1.5
//realloc sparse matrix NZ arrays
inline int reallocSpMatrix(spmat* mat,uint newSize){
    mat->NZ *= newSize;
    void* tmp;
    if (!(tmp = realloc(mat->AS,mat->NZ * sizeof(*(mat->AS))))){
        ERRPRINT("reallocSpMatrix:  realloc AS errd\n");
        return EXIT_FAILURE;
    }
    mat->AS = tmp;
    if (!(tmp = realloc(mat->JA,mat->NZ * sizeof(*(mat->JA))))){
        ERRPRINT("reallocSpMatrix:  realloc JA errd\n");
        return EXIT_FAILURE;
    }
    mat->JA = tmp;
    return EXIT_SUCCESS;
}
*/
////MISC
//print useful information about 3SPGEMM about to compute
void print3SPGEMMCore(spmat* R,spmat* AC,spmat* P,CONFIG* conf);
#endif
