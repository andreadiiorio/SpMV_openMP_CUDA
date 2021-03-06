/*
Copyright Andrea Di Iorio 2022
This file is part of SpMV_OMP_CUDA
SpMV_OMP_CUDA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SpMV_OMP_CUDA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SpMV_OMP_CUDA.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include "sparseMatrix.h"
#include "macros.h"
#include "utils.h"

///inline export here 
//SPMV_CHUNKS_DISTR spmvChunksFair; 
spmat* allocSpMatrix(ulong rows, ulong cols);
int allocSpMatrixInternal(ulong rows, ulong cols, spmat* mat);
void freeSpmatInternal(spmat* mat);
void freeSpmat(spmat* mat);

////////////////////////  CSR SPECIFIC -- TODO RENAME //////////////////
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
        offsets[ IDX2D(r,0,gridCols) ] = j;  //row's first gc start is costrained
        //navigate column groups inside current row
        for (ulong gc=1,gcStartCol;  gc<gridCols;  gc++){
            gcStartCol = UNIF_REMINDER_DISTRI_STARTIDX(gc,_colBlock,_colBlockRem);
            //goto GroupCols start entry,keeping A's nnz entries navigation (idx j)
            //for (ulong c=A->JA[j]; c<gcStartCol && j < A->IRP[r+1]; c=A->JA[++j]);
            while ( j < A->IRP[r+1] &&  A->JA[j] < gcStartCol )  j++;
            offsets[ IDX2D(r,gc,gridCols) ] = j;  //row's gc group startIdx
        }
    }
    offsets[subRowsN] = A->NZ;  //last row's partition end costrained
    return offsets;
}

spmat* colsPartitioningUnifRanges(spmat* A,ulong gridCols){
    spmat *colParts, *colPart;
    ulong _colBlock = A->N/gridCols, _colBlockRem = A->N%gridCols, *colPartsLens=NULL, *tmpJA;
    double* tmpAS;
    ///alloc/init partitions structures
    if (!(colParts = calloc(gridCols, sizeof(*colParts)))){
        ERRPRINT("colsPartitioningUnifRanges\tcolumns partitions of A calloc fail\n");
        return NULL;
    }
    for (ulong i=0,colBlock; i<gridCols; i++){
        colBlock = UNIF_REMINDER_DISTRI(i,_colBlock,_colBlockRem);
        colPart  = colParts + i;
        if (allocSpMatrixInternal(A->M,colBlock,colPart)){
            ERRPRINT("colsPartitioningUnifRanges\tallocSpMatrixInternal partition err\n");
            goto _err;
        }
        //TODO TODO overalloc A cols partitions NZ arrays, then realloc
        if (!(colPart->AS = malloc(A->NZ * sizeof(*A->AS)))){
            ERRPRINT("colPart of A overalloc of AS errd\n");
            goto _err;
        }
        if (!(colPart->JA = malloc(A->NZ * sizeof(*A->JA)))){
            ERRPRINT("colPart of A overalloc of JA errd\n");
            goto _err;
        }
    }
    //for each A col partition -> last copied nz index = nnz copied ammount
    colPartsLens = calloc(gridCols, sizeof(colPartsLens));
    if (!colPartsLens){
        ERRPRINT("colsPartitioningUnifRanges: colPartsLens calloc errd\n");
        goto _err;
    }
    //OFFSET BASED COPY OF A.COL_GROUPS -> O( A.NZ )
    /* TODO
     * Parallelize: 2for collapse OMP, gcEndCol -> startIdxize, ...
     * oppure wrappare cio in static inline 
     * ed implementazione alternativa: implementazione basata su offset da colsOffsetsPartitioningUnifRanges  
     */
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
    free(colPartsLens);
    return colParts;
    _err:
    for (ulong i=0; i<gridCols; i++)   freeSpmatInternal(colParts+i);
    if(colParts)        free(colParts);
    if(colPartsLens)    free(colPartsLens);
    return NULL;
}

////////////////////////	ELL AUX FUNCS //////////////////////////////
spmat* ellTranspose(spmat* m){
	//alloc components
	spmat* out = calloc(1,sizeof(*out));
	if (!out){
		ERRPRINT("ellTranspose: out callc errd\n");
		return NULL;
	}
	if (!(out->AS = malloc(sizeof(*(out->AS)) * m->M*m->MAX_ROW_NZ))){
		ERRPRINT("ellTranspose: invalid malloc AS\n");
		goto _err;
	}
	if (!(out->JA = malloc(sizeof(*(out->JA)) * m->M*m->MAX_ROW_NZ))){
		ERRPRINT("ellTranspose: invalid malloc JA\n");
		goto _err;
	}
	#ifdef ROWLENS
	if (!(out->RL = malloc(sizeof(*(out->RL)) * m->M))){
		ERRPRINT("ellTranspose: invalid malloc RL\n");
		goto _err;
	}
	memcpy(out->RL,m->RL,sizeof(*(out->RL)) * m->M);
	#endif
	//init
	out -> NZ 			= m->NZ;
	out -> M  			= m->MAX_ROW_NZ;
	out -> N  			= m->M;
	out -> MAX_ROW_NZ 	= m->M;
	//transposed copy
	for(ulong r=0; r<m->M; r++){
		for(ulong c=0; c<m->MAX_ROW_NZ; c++){
			out->JA[IDX2D(c,r,out->MAX_ROW_NZ)] = m->JA[IDX2D(r,c,m->MAX_ROW_NZ)];
			out->AS[IDX2D(c,r,out->MAX_ROW_NZ)] = m->AS[IDX2D(r,c,m->MAX_ROW_NZ)];
		}
	}
	VERBOSE	printf("Transposed ELL mat at %p in %p\n",m,out);
	return out;
	
	_err:
	freeSpmat(out);
	return NULL;
}
////////////////////////////////////////////////////////////////////////
int spmatDiff(spmat* A, spmat* B){
    if (A->NZ != B->NZ){
        ERRPRINT("NZ differ\n");
        return EXIT_FAILURE;
    }
    if (doubleVectorsDiff(A->AS,B->AS,A->NZ,NULL)){
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


///unit test embbeded
#ifdef SPARSEUTILS_MAIN_TEST

////INTERNAL TEST FUNCTIONS
//test that each row's partition from colsOffsetsPartitioningUnifRanges is in the correct index range
#include <alloca.h>
int testColsOffsetsPartitioningUnifRanges(spmat* mat,ulong gridCols,ulong* partsOffs){
    ulong _colBlock = mat->N/gridCols, _colBlockRem = mat->N%gridCols;
    ulong j=0;    //CSR scanning nnz idx
    ulong* colPartsPopulations = alloca(gridCols * sizeof(*colPartsPopulations));
    memset(colPartsPopulations,0,gridCols * sizeof(*colPartsPopulations));
    for (ulong r=0,pId=0; r<mat->M; r++){
        for (ulong gc=0,pStartIdx,pEndIdx; gc<gridCols; gc++,pId++){
            pStartIdx = UNIF_REMINDER_DISTRI_STARTIDX(gc,_colBlock,_colBlockRem);
            pEndIdx   = UNIF_REMINDER_DISTRI_STARTIDX(gc+1,_colBlock,_colBlockRem)-1; 
            //pId=IDX2D(r,gc,gridCols);
            for (ulong idx=partsOffs[pId],c; idx<partsOffs[pId+1]; idx++,j++){
                c = mat->JA[idx];
                assert( j == idx ); //consecutive index in partitioning
                assert( pStartIdx <= c && c <= pEndIdx );               //colRange
                assert( mat->IRP[r] <= idx && idx <= mat->IRP[r+1] );   //rowRange
            }
            colPartsPopulations[gc] += partsOffs[pId+1] - partsOffs[pId]; 
        }
    }
    assert( j == mat->NZ );
    ulong s=0;
    for (ulong gc=0,partSize; gc < gridCols; gc++,s+=partSize){
        partSize = colPartsPopulations[gc];
        double partShare=partSize/(double)mat->NZ,partsAvg=1/(double)gridCols;
        double partShareAvgDiff = partShare - partsAvg;
        printf("colPartition %lu has:\t%lu = %lf of NNZ\t\t .. %lf\tAVG diff\n",
          gc,partSize,partShare,partShareAvgDiff);
    }
    assert( s == mat->NZ ); //TODO DUPLICATED
    return EXIT_SUCCESS;
}

CONFIG Conf = {
    .gridRows = 8,
    .gridCols = 8,
};

#include "parser.h"
int main(int argc, char** argv){
    int out=EXIT_FAILURE;
    if (init_urndfd())  return out;
    if (argc < 2 )  {ERRPRINT("COO MATRIX FOR TEST"); return out;}
    ////parse sparse matrix and dense vector
    spmat* mat;
    char* trgtMatrix = TMP_EXTRACTED_MARTIX;
    if (extractInTmpFS(argv[1],TMP_EXTRACTED_MARTIX) < 0)   trgtMatrix = argv[1];
    if (!(mat = MMtoCSR(trgtMatrix))){
        ERRPRINT("err during conversion MM -> CSR\n");
        return out;
    }
    ////partitioning test
    ulong* colsPartitions = colsOffsetsPartitioningUnifRanges(mat,Conf.gridCols);
    if (!colsPartitions)    goto _free;
    if (testColsOffsetsPartitioningUnifRanges(mat,Conf.gridCols,colsPartitions))  goto _free;

    out=EXIT_SUCCESS;
    printf("testColsOffsetsPartitioningUnifRanges passed with "
           "mat: %lux%lu-%luNNZ\tgrid: %dx%d\n",
            mat->M,mat->N,mat->NZ,Conf.gridRows,Conf.gridCols);
    _free:
    if (colsPartitions) free(colsPartitions);

    return out;
}
#endif
