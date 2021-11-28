#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sparseMatrix.h"
#include "mmio.h"
#include "parser.h"
#include "macros.h"
#include "utils.h"


double* readVector(char* fpath,ulong* size){
    double *out,*tmp;
    ulong i=0,vectorSize=VECTOR_STEP_MALLOC;
    FILE* fp = fopen(fpath,"r");
    if (!fp){
        perror("fopen vector file");
        return NULL;
    }
    if (!(out = malloc(VECTOR_STEP_MALLOC*sizeof(*out)))){ 
        ERRPRINT("vector read malloc fail for file\n");
        return NULL;
    }
    while (1){
        if (i >= vectorSize ){ //realloc the array
            if (!(tmp=realloc(out,vectorSize+VECTOR_STEP_MALLOC*sizeof(*out)))){
                ERRPRINT("realloc errd\n");
                goto _err;
            }
            out = tmp;
        }
        if (fscanf(fp,"%lf ",out + i++ ) != 1){
            perror("invalid fscanf");
            goto _err;
        }
    }
    //REALLOC THE ARRAY TO THE FINAL SIZE
    *size = --i;
    if (!(tmp = realloc(out,*size*sizeof(*out)))){
        ERRPRINT("realloc errd\n");
        goto _err;
    }
    VERBOSE printf("readed array rellocd to %lu bytes\n",*size);
    out = tmp;
    
    return tmp;
    _err:
    free(out);
    return NULL;
}
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

/* parse MatrixMarket matrix in @fp into @mat
 * expand simmetric matrix into a normal matrix with both parts
 * return sparse matrix
 */
static int _MMtoCSR(spmat* mat, FILE *fp, MM_typecode mcode){
    int out=EXIT_FAILURE, scanndRet=0;
    ulong nzIdx=0;         //expanded num of nz (in case of sym matrix)
    ulong idx;
    double val;         //current entry val
    ulong row,col;       //current entry's row,col from MM -> 1 based
    ulong* rowsNextCol = NULL;  //for each row -> next entry progressive idx
    int* _rowsLastCol = NULL; //for each row -> last added entry's columnIdx 
    entry* entries = NULL;     //MM parsed and exmpanded entries
    ///consistency checks
    ///init
    if (mm_is_symmetric(mcode)) {
        mat->NZ *= 2;
    }
    if (!(entries         = malloc(mat->NZ * sizeof(*entries))))  goto _free;
        if (!(rowsNextCol = calloc(mat->M,sizeof(*rowsNextCol)))) goto _free;
    CONSISTENCY_CHECKS{
        if (!(_rowsLastCol = malloc(mat->M*sizeof(*rowsNextCol)))) goto _free;
        memset(_rowsLastCol,-1,mat->M*sizeof(*_rowsLastCol));
    }
    ///parse MM fp lines
    while (1) { // Reading the fp until EOF
        if          (mm_is_pattern(mcode)){
            scanndRet = fscanf(fp, "%lu %lu\n", &row, &col);
            val = 1.0;
        } else if   (mm_is_real(mcode) || (mm_is_integer(mcode))){
            scanndRet = fscanf(fp, "%lu %lu %lf\n", &row, &col, &val);
        }
        

        if (scanndRet == EOF){ //TODO more strict check with type&ret?
            if (ferror(fp)){
                perror("fscanf EOF");
                goto _free;
            } else  break;
        }
        CONSISTENCY_CHECKS{ 
            //NNZERO VAL CHECKS 
            //if (!val){    //TODO SOMETIMES THERE IS IN SOME MATRIXES...
            //    ERRPRINT("invalid sparse matrix with 0 entry explicitly stored\n");
            //    goto _free;
            //}
            //TODO USELESS ? ? ?
            if ((mm_is_pattern(mcode) && scanndRet != 2) || 
                (!mm_is_pattern(mcode) && scanndRet != 3)){
                ERRPRINT("invalid matrix: not consistent entry scannable\n");
                goto _free;
            }
        }
       ////ADD THE CURRENT MATRIX ENTRY
       //IRP will be builded as definition of Row start Pointers
       //from shifted row lenghts: in IRP[1] <- len Row[0], IRP[2] <- len ROw[1]
       //after, cumulating the row lens as increment to get row start idx ptrs
       mat -> IRP[row]++;   //now just count entries per row 1based
       entries[nzIdx++]=(entry) { .row=row-1, .col=col-1, .val=val };

        //also mirrored entry if sym.matrix with reflected idx inside matrix limits
        if (mm_is_symmetric(mcode) && 
            row != col && row <= mat->N && col <= mat->M ){

            swap(row,col);
            mat-> IRP[row]++;   //now just count entries per row
            entries[nzIdx++]=(entry) { .row=row-1, .col=col-1, .val=val };
        }
    }
    /*if (SORT_MM_ENTRIES){ //force entries sorting before pack in CSR 
        //TODO SORT entries by key=row,col 
    }*/
    
    ////build core struct of CSR
    mat->NZ = nzIdx; //all rappresented element (expanding symmetric ones)
    
    if(!(mat->JA = malloc(nzIdx*sizeof(*(mat->JA))))){
        ERRPRINT("JA malloc err\n");
        goto _free;
    }
    if(!(mat->AS = malloc(nzIdx*sizeof(*(mat->AS))))){
        ERRPRINT("AS malloc err\n");
        goto _free;
    }
#ifdef ROWLENS
    //auxiliary row lengths array copied from original IRP just after parse
    //exploit shifted row lenghts previously computed
    memcpy(mat->RL,mat->IRP + 1,sizeof(*mat->IRP) * mat->M);
    //for (ulong i=1;i<mat->M+1;i++)mat->RL[i-1] = mat->IRP[i]-mat->IRP[i-1];//TODO usable with definition IRP
#endif
    //IRP: trasform rows lens as increments to build row index "pointer"
    //0th -> 0 mandatory; 1th = 0th row len, ...., M+1th = end of Mth row
    for (ulong i=2; i<mat->M+1; i++)    mat->IRP[i] += mat->IRP[i-1];
    CONSISTENCY_CHECKS{
        if (mat->IRP[mat->M] != mat->NZ){
            fprintf(stderr,"IRP[M] %lu != NZ %lu\n",mat->IRP[mat->M],mat->NZ);
            goto _free;
        }
    }
    for (ulong i=0; i<nzIdx; i++) {
        row = entries[i].row;
        col = entries[i].col;
        val = entries[i].val;
        CONSISTENCY_CHECKS{
            if (_rowsLastCol[row] >= (int) col){
                ERRPRINT("not sorted rows\n");
                goto _free;
            }
            _rowsLastCol[row] = col;
        }
        idx = mat -> IRP[row] + rowsNextCol[row]++;
        mat -> AS[idx] = val;
        mat -> JA[idx] = col;
    }
    out = EXIT_SUCCESS;
    //goto _free;

    _free:
    if (entries)        free(entries);
    if (rowsNextCol)    free(rowsNextCol);
    CONSISTENCY_CHECKS{ if(_rowsLastCol)    free(_rowsLastCol);}

    return out;
}

spmat* MMtoCSR(char* matPath){
    MM_typecode mcode;
    spmat* mat = NULL;
    FILE* fp = fopen(matPath, "r");
    if (!fp){
        perror("fopen");
        goto err;
    }
    //banner -> parse  matrix specs
    if (mm_read_banner(fp, &mcode) != 0) {
        fprintf(stderr,"mm_read_banner err at:%s\n",matPath);
        goto err;
    }
    //assert matrix is compatible with this app scope
    if (MMCheck(mcode))     goto out;
    //alloc sparse matrix components
    if (!(mat = calloc(1,sizeof(*mat)))){
        ERRPRINT(" mat struct alloc errd");
        goto err;
    }
    //parse sizes
    //TODO OVERCOME uint limitation?
    if(mm_read_mtx_crd_size(fp, (uint*) &mat->M, (uint*) &mat->N, (uint*) &mat->NZ)){
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

    if (_MMtoCSR(mat, fp, mcode)){
        fprintf(stderr,"MAT PARSE TO CSR ERR at:%s\n",matPath);
        goto err;
    }
    
    out:
    fclose(fp);
    return mat;

    err:
    if (mat)    freeSpmat(mat);
    if (fp)         fclose(fp);
    return NULL;
}
//TODO spmat* MMtoHELL(char* matPath){}
