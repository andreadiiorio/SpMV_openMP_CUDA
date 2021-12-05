#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sparseMatrix.h"
#include "mmio.h"
#include "parser.h"
#include "macros.h"
#include "utils.h"

int MMCheck(MM_typecode mcode) {
    if (!mm_is_matrix(mcode)){  //consistency checks among flags in @mcode
        ERRPRINT("invalid matrix: not a matrix\n");
        return EXIT_FAILURE;
    }
    if (mm_is_dense(mcode) ){   //|| mm_is_array(mcode) ){
        ERRPRINT("invalid matrix: not a supported sparse matrix\tDENSE MAT\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

entry* MMtoCOO(ulong* NZ, FILE *fp, MM_typecode mcode){
    int scanndRet=0;
    ulong nzIdx=0;                //expanded num of nz (in case of sym matrix)
    ulong row,col;                //current entry's row,col from MM -> 1 based
    double val;
    entry* entries = NULL;        //COO parsed entries
    ///init
    if (mm_is_symmetric(mcode)){
        (*NZ) *= 2;
        VERBOSE     printf("simmetric matrix\n");
    }
    if (!(entries     = malloc(*NZ * sizeof(*entries)))){
        ERRPRINT("MMtoCOO:  entries malloc errd\n");
        return NULL;
    }
    ///parse MM fp lines into COOordinate entries
    while (1) { // Reading the fp until EOF
        if (mm_is_pattern(mcode)){
            scanndRet = fscanf(fp, "%lu %lu\n", &row, &col);
            val = 1.0;
        } else if   (mm_is_real(mcode) || (mm_is_integer(mcode))){
            scanndRet = fscanf(fp, "%lu %lu %lf\n", &row, &col, &val);
        }
        

        if (scanndRet == EOF){ //TODO more strict check with type&ret?
            if (ferror(fp)){
                perror("fscanf EOF");
                goto _err; 
            } else  break;
        }
        CONSISTENCY_CHECKS{ 
            //TODO USELESS ? ? ?
            if ((mm_is_pattern(mcode) && scanndRet != 2) || 
                (!mm_is_pattern(mcode) && scanndRet != 3)){
                ERRPRINT("invalid matrix: not consistent entry scannable\n");
                goto _err; 
            }
        }
       ////ADD THE CURRENT MATRIX ENTRY
       //TODO OLD VERSION IRP will be builded as definition of Row start Pointers
       //from shifted row lenghts: in IRP[1] <- len Row[0], IRP[2] <- len ROw[1]
       //after, cumulating the row lens as increment to get row start idx ptrs
       //mat -> IRP[row]++;   //now just count entries per row 1based
       entries[nzIdx++]=(entry) { .row=row-1, .col=col-1, .val=val };

        //also mirrored entry if sym.matrix with reflected idx inside matrix limits
        if (mm_is_symmetric(mcode) && row != col ){
            //TODO COSTRAINED FORMAT ?&& row <= mat->N && col <= mat->M ){

            swap(row,col);
            //mat-> IRP[row]++;   //now just count entries per row
            entries[nzIdx++]=(entry) { .row=row-1, .col=col-1, .val=val };
        }
    }
    /*if (SORT_MM_ENTRIES){ //force entries sorting before pack in CSR 
        //TODO SORT entries by key=row,col 
    }*/
   
    CONSISTENCY_CHECKS  assert( nzIdx == *NZ );
    return entries;
    
    _err:
    free(entries);  return NULL;
} 

int COOtoCSR(entry* entries, spmat* mat){ 
    int out = EXIT_FAILURE;
    ulong idx;
    long* _rowsLastCol = NULL;    //for each row -> last added entry's columnIdx 
    ulong* rowsNextIdx = NULL;    //for each row -> next entry progressive idx
    if (!(rowsNextIdx = calloc(mat->M,sizeof(*rowsNextIdx)))){
        ERRPRINT("MMtoCOO:  rowsNextIdx calloc errd\n");
        goto _end;
    }
    CONSISTENCY_CHECKS{ //alloc and init aux arr for entries sort check
        if (!(_rowsLastCol = malloc(mat->M*sizeof(*_rowsLastCol)))){
            ERRPRINT("MMtoCOO:  _rowsLastCol malloc errd\n");
            goto _end;
        }
        memset(_rowsLastCol,-1,mat->M*sizeof(*_rowsLastCol));
    }
    //get rowLens -> IRP (partial)
    for (ulong i=0; i<mat->NZ; i++)     mat->IRP[entries[i].row+1]++;
    #ifdef ROWLENS
    memcpy(mat->RL,mat->IRP + 1,sizeof(*mat->IRP) * mat->M);
    #endif
    //IRP: trasform rows lens as increments to build row index "pointer"
    //0th -> 0 mandatory; 1th = 0th row len, ...., M+1th = end of Mth row
    for (ulong i=2; i<mat->M+1; i++)    mat->IRP[i] += mat->IRP[i-1];
    CONSISTENCY_CHECKS  assert(mat->IRP[mat->M] == mat->NZ);
    //entries write in CSR format
    entry* e;
    for (ulong i=0; i<mat->NZ; i++) {
        e = entries+i;
        CONSISTENCY_CHECKS{ //TODO CHECK IF COO ENTRIES ARE SORTED
            #pragma message("COO sorting check enabled")
            if (_rowsLastCol[e->row] >= (long) e->col){
                ERRPRINTS("not sorted entry:%ld,%ld,%lf",e->row,e->col,e->val);
                goto _end;
            }
            _rowsLastCol[e->row] = e->col;
        }
        idx = mat -> IRP[e->row] + rowsNextIdx[e->row]++;
        mat -> AS[idx] = e->val;
        mat -> JA[idx] = e->col;
    }
    
    out = EXIT_SUCCESS;

    _end:
    if(rowsNextIdx)                         free(rowsNextIdx);
    CONSISTENCY_CHECKS{ if(_rowsLastCol)    free(_rowsLastCol);}

    return out;
}

spmat* MMtoCSR(char* matPath){
    MM_typecode mcode;
    spmat* mat = NULL;
    entry* entries=NULL;
    FILE* fp = fopen(matPath, "r");
    if (!fp){
        perror("fopen");
        return NULL;
    }
    //banner -> parse  matrix specs
    if (mm_read_banner(fp, &mcode) != 0) {
        fprintf(stderr,"mm_read_banner err at:%s\n",matPath);
        goto err;
    }
    //assert matrix is compatible with this app scope
    if (MMCheck(mcode))     goto err;
    //alloc sparse matrix components
    if (!(mat = calloc(1,sizeof(*mat)))){
        ERRPRINT(" mat struct alloc errd");
        goto err;
    }
    //parse sizes
    //TODO OVERCOME uint limitation?
    if(mm_read_mtx_crd_size(fp,(uint*) &mat->M, (uint*) &mat->N, (uint*) &mat->NZ)){
        fprintf(stderr,"mm_read_mtx_crd_size err at %s:\n",matPath);
        goto err;
    }
    if (!(mat->IRP = calloc(mat->M+1,sizeof(*(mat->IRP))))){
        ERRPRINT("IRP calloc err\n");
        goto err;
    }
#ifdef ROWLENS
    if (!(mat->RL = calloc(mat->M,sizeof(*(mat->RL))))){
        ERRPRINT("IRP calloc err\n");
        goto err;
    }
#endif

    if (!(entries = MMtoCOO(&(mat->NZ), fp, mcode)) ){
        ERRPRINTS("MAT PARSE TO CSR ERR at:%s\n",matPath);
        goto err;
    }
    ////alloc core struct of CSR
    if(!(mat->JA = malloc(mat->NZ*sizeof(*(mat->JA))))){
        ERRPRINT("JA malloc err\n");
        goto err; 
    }
    if(!(mat->AS = malloc(mat->NZ*sizeof(*(mat->AS))))){
        ERRPRINT("AS malloc err\n");
        goto err;  
    }
    if (COOtoCSR(entries,mat))  goto err;
   
    goto _free;


    err:
    if (mat)        freeSpmat(mat);
    mat = NULL;
    _free:
    fclose(fp);
    if (entries)    free(entries);
    return mat;
}
//TODO spmat* MMtoHELL(char* matPath){}
