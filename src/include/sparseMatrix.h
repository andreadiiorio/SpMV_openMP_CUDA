//sparse matrix def & aux
//TODO adapt to work on both CUDA->ELL and std CSR
#ifndef SPARSEMATRIX
#define SPARSEMATRIX 

#include "macros.h"
#include "config.h"
typedef struct{
    ulong NZ,M,N;
    ulong* JA;
    //CSR SPECIFIC
#ifdef ROWLENS
    ulong* RL;   //row lengths
#endif
    ulong* IRP;
	//ELL SPECIFIC
    ulong MAX_ROW_NZ;

    double *AS; 
	#ifdef __CUDACC__	//CUDA SPECIFIC
	//pitch of ELL* mat allocated as number of element of respective arrays
	size_t pitchJA;
	size_t pitchAS;
	#endif
} spmat; //describe a sparse matrix

////Sparse vector accumulator -- corresponding to a matrix portion
typedef struct{
    //ulong    r;     //row index in the corresponding matrix
    //ulong    c;     //col index in the corresponding matrix
    ulong   len;   //rowLen
    double* AS;    //row nnz    values
    ulong*  JA;    //row nnz    colIndexes
} SPACC; 

//TODO ASSERT LEN>0 ommitted
inline int BISECT_ARRAY(ulong target, ulong* arr, ulong len){
    //if (len == 0)              return FALSE;
    if (len <= 1)              return *arr == target; 
    ulong middleIdx = (len-1) / 2;  //len=5-->2, len=4-->1
    ulong middle    = arr[ middleIdx ];
    if      (target == middle)  return TRUE;
    else if (target <  middle)  return BISECT_ARRAY(target,arr,middleIdx); 
    else    return BISECT_ARRAY(target,arr+middleIdx+1,middleIdx + (len-1)%2);
}

/*
 * return !0 if col @j idx is in row @i of sparse mat @smat
 * bisection used --> O(log_2(ROWLENGHT))
 */
inline int IS_NNZ(spmat* smat,ulong i,ulong j){
    ulong rStart = smat->IRP[i];
    ulong rLen   = smat->IRP[i+1] - rStart;
    if (!rLen)  return FALSE;
    return BISECT_ARRAY(j,smat->JA + rStart,rLen);
}
inline int IS_NNZ_linear(spmat* smat,ulong i,ulong j){    //linear -> O(ROWLENGHT)
    int out = 0;
    for (ulong x=smat->IRP[i]; x<smat->IRP[i+1] && !out; x++){
        out = (j == smat->JA[x]); 
    } 
    return out;
}
////aux functions
//free sparse matrix
inline void freeSpmatInternal(spmat* mat){
    free(mat->AS);  
    free(mat->JA);  
    free(mat->IRP);  
#ifdef ROWLENS
    free(mat->RL);
#endif 
}
inline void freeSpmat(spmat* mat){
    freeSpmatInternal(mat);
    free(mat);
}

//free max aux structs not NULL pointed
inline void freeSpAcc(SPACC* r){ 
    free(r->AS);
    free(r->JA);
}
////alloc&init functions
//alloc&init internal structures only dependent of dimensions @rows,@cols
inline int allocSpMatrixInternal(ulong rows, ulong cols, spmat* mat){
    mat -> M = rows;
    mat -> N = cols;
    if (!(mat->IRP=(typeof(mat->IRP)) calloc(mat->M+1,sizeof(*(mat->IRP))))){ //0set only for 0th
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
inline spmat* allocSpMatrix(ulong rows, ulong cols){

    spmat* mat;
    if (!(mat = (typeof(mat)) calloc(1,sizeof(*mat)))) { 
        ERRPRINT("mat  calloc failed\n");
        return NULL;
    }
    if (allocSpMatrixInternal(rows,cols,mat)){
        free(mat);
        return NULL;
    }
    return mat;
}

//////////////////////////////// CSR SPECIFIC /////////////////////////////////
///SPARSE MATRIX PARTITIONING
/*
 * partition CSR sparse matrix @A in @gridCols columns partitions 
 * returning an offsets matrix out[i][j] = start of jth colPartition of row i
 * subdivide @A columns in uniform cols ranges in the output 
 */
ulong* colsOffsetsPartitioningUnifRanges(spmat* A,ulong gridCols);

/*
 * partition CSR sparse matrix @A in @gridCols columns partitions as 
 * indipended and allocated sparse matrixes and return them
 * subdivide @A columns in uniform cols ranges in the output 
 */
spmat* colsPartitioningUnifRanges(spmat* A,ulong gridCols);
////////////////////////	ELL AUX FUNCS //////////////////////////////
spmat* ellTranspose(spmat* m);
///////////////////////////////////////////////////////////////////////////////

/*  
    check if sparse matrixes A<->B differ up to 
    DOUBLE_DIFF_THREASH per element
*/
int spmatDiff(spmat* A, spmat* B);
////dyn alloc of spGEMM output matrix
/*
///size prediction of AB = @A * @B
inline ulong SpGEMMPreAlloc(spmat* A,spmat* B){
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
inline int reallocSpMatrix(spmat* mat,ulong newSize){
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
#endif
