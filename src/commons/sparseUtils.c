#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include "macros.h"
#include "sparseMatrix.h"
#include "utils.h"

///SPARSE MATRIX PARTITIONING
ulong* colsOffsetsPartitioningUnifRanges(spmat* A,ulong gridCols){
    ulong subRowsN = A->M * gridCols;
    ulong _colBlock = A->N/gridCols, _colBlockRem = A->N%gridCols;
    ulong* offsets = malloc( (subRowsN+1) * sizeof(*offsets) );
    if (!offsets)  {
        ERRPRINT("colsOffsetsPartitioningUnifRanges:\toffsets malloc errd\n");
        return NULL;
    }
    ///OFFSETS COMPUTE FOR COL GROUPS -> O( A.NZ )
    for (ulong r=0, j=0;     r<A->M;     j=A->IRP[++r]){
        //navigate column groups inside current row
        for (ulong gc=0,gcStartCol=0;  gc<gridCols;  gc++){
            //goto GroupCols start entry,keeping A's nnz entries navigation (idx j)
            //for (ulong c=A->JA[j]; c<gcStartCol && j < A->IRP[r+1]; c=A->JA[++j]);
            while ( j < A->IRP[r+1] &&  A->JA[j] < gcStartCol )  j++;
            offsets[ IDX2D(r,gc,gridCols) ] = j;  //row's gc group startIdx
            gcStartCol += UNIF_REMINDER_DISTRI(gc,_colBlock,_colBlockRem);
        }
    }
    offsets[subRowsN] = A->NZ;
    return offsets;
}

spmat* colsPartitioningUnifRanges(spmat* A,ulong gridCols){
    spmat *colParts, *colPart;
    ulong _colBlock = A->N/gridCols, _colBlockRem = A->N%gridCols, *tmpJA;
    double* tmpAS;
    ///alloc/init partitions structures
    if (!(colParts = calloc(gridCols, sizeof(*colParts)))){
        ERRPRINT("colsPartitioningUnifRanges\tcolumns partitions of A calloc fail\n");
        return NULL;
    }
    for (ulong i=0,colBlock; i<gridCols; i++){
        colBlock = UNIF_REMINDER_DISTRI(i,_colBlock,_colBlockRem);
        if (allocSpMatrixInternal(A->M,colBlock,colParts + i)){
            ERRPRINT("colsPartitioningUnifRanges\tallocSpMatrixInternal partition err\n");
            goto _err;
        }
        //TODO overalloc A cols partitions NZ arrays, then realloc
        if (!((colParts+i)->AS = malloc(A->NZ * sizeof(*A->AS)))){
            ERRPRINT("colPart of A overalloc of AS errd\n");
            goto _err;
        }
        if (!((colParts+i)->JA = malloc(A->NZ * sizeof(*A->JA)))){
            ERRPRINT("colPart of A overalloc of JA errd\n");
            goto _err;
        }
    }
    //for each A col partition -> last copied nz index = nnz copied ammount
    ulong* colPartsLens = alloca(gridCols * sizeof(colPartsLens));
    memset(colPartsLens, 0, sizeof(*colPartsLens) * gridCols);
    //OFFSET BASED COPY OF A.COL_GROUPS -> O( A.NZ )
    for (ulong r=0, j=0;     r<A->M;     j=A->IRP[++r]){
        //navigate column groups inside current row
        for (ulong gc=0,gcEndCol=0,i;  gc<gridCols ;  gc++,j+=i){
            i = 0;  //@i=len current subpartition of row @r to copy
            colPart = colParts + gc;
            colPart->IRP[r] = colPartsLens[gc];
            gcEndCol += UNIF_REMINDER_DISTRI(gc,_colBlock,_colBlockRem);
            //goto next GroupCols,keeping A's nnz entries navigation ( index j+i )
            //for (ulong c=A->JA[j+i]; c<gcEndCol && j+i  < A->IRP[r+1]; c=A->JA[j+ ++i]);
            while ( j+i < A->IRP[r+1] && A->JA[j+i] < gcEndCol ) i++;
            memcpy(colPart->AS+colPart->IRP[r], A->AS+j, i*sizeof(*(A->AS)));
            memcpy(colPart->JA+colPart->IRP[r], A->JA+j, i*sizeof(*(A->JA)));
            
            colPartsLens[gc] += i;
#ifdef ROWLENS
            colPart->RL[r] = i;
#endif
        }
    }
    //TODO realloc overallcd A parts NZ arrays (TODO -> downsizing -> nofails?)
    for (ulong i=0,partLen; i<gridCols; i++){
        colPart = colParts + i;
        partLen = colPartsLens[i];
        if (!(tmpAS = realloc(colPart->AS,partLen*sizeof(*(colPart->AS))))){
            ERRPRINT("realloc overallocated cols partition AS array\n");
            goto _err;
        }
        colPart->AS = tmpAS;
        if (!(tmpJA = realloc(colPart->JA,partLen*sizeof(*(colPart->JA))))){
            ERRPRINT("realloc overallocated cols partition JA array\n");
            goto _err;
        }
        colPart->JA         = tmpJA;
        colPart->NZ         = partLen;
        colPart->IRP[A->M]  = partLen;
    }
    return colParts;
    _err:
    for (ulong i=0; i<gridCols; i++)   freeSpmatInternal(colParts+i);
    return NULL;
}

int spmatDiff(spmat* A, spmat* B){
    if (A->NZ != B->NZ){
        ERRPRINT("NZ differ\n");
        return EXIT_FAILURE;
    }
    if (doubleVectorsDiff(A->AS,B->AS,A->NZ)){
        ERRPRINT("AS DIFFER\n");
        return EXIT_FAILURE;
    }
    if (memcmp(A->JA,B->JA,A->NZ)){
        ERRPRINT("JA differ\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

double* CSRToDense(spmat* sparseMat){
    double* denseMat;
    ulong i,j,idxNZ,denseSize;
    if (__builtin_umull_overflow(sparseMat->M,sparseMat->N,&denseSize)){
        ERRPRINT("overflow in dense allocation\n");
        return NULL;
    }
    if (!(denseMat = calloc(denseSize, sizeof(*denseMat)))){
        ERRPRINT("dense matrix alloc failed\n");
        return  NULL;
    }
    for (i=0;i<sparseMat->M;i++){
        for (idxNZ=sparseMat->IRP[i]; idxNZ<sparseMat->IRP[i+1]; ++idxNZ){
             j = sparseMat->JA[idxNZ];
             //converting sparse item into dense entry
             denseMat[(ulong) IDX2D(i,j,sparseMat->N)] = sparseMat->AS[idxNZ]; 
        }
    }
    return denseMat;
}
void printSparseMatrix(spmat* spMatrix,char justNZMarkers){
    double* denseMat = CSRToDense(spMatrix);
    if (!denseMat)  return;
    printMatrix(denseMat,spMatrix->M,spMatrix->N,justNZMarkers);
    free(denseMat);
}

void print3SPGEMMCore(spmat* R,spmat* AC,spmat* P,CONFIG* conf){
    printf("@COARSENING AC: %lux%lu ---> %lux%lu\tconf grid: %ux%u,\tNNZ:%lu-%lu-%lu\t",
      AC->M,AC->N, R->M,P->N, conf->gridRows,conf->gridCols, R->NZ,AC->NZ,P->NZ);
}
