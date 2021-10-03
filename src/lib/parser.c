#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sparseMatrix.h"
#include "mmio.h"
#include "parser.h"
#include "utils.h"


double* readVector(char* fpath){
    double *out,*tmp;
    uint i=0,vectorSize=VECTOR_STEP_MALLOC;
    FILE* fp = fopen(fpath,"r");
    if (!fp){
        perror("fopen vector file");
        return NULL;
    }
    if (!(out = malloc(VECTOR_STEP_MALLOC*sizeof(*out)))){ 
        fprintf(stderr,"vector read malloc fail for file\n");
        return NULL;
    }
    while (1){
        if (i >= vectorSize ){ //realloc the array
            if (!(tmp = reallocarray(out,
                vectorSize+VECTOR_STEP_MALLOC,sizeof(*out)))){
                fprintf(stderr,"reallocarray errd\n");
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
    if (!(tmp = reallocarray(out,--i,sizeof(*out)))){
        fprintf(stderr,"reallocarray errd\n");
        goto _err;
    }
    out = tmp;
    
    return tmp;
    _err:
    free(out);
    return NULL;
}
int MMCheck(MM_typecode mcode) {
    if (!mm_is_matrix(mcode)){
        fprintf(stderr,"invalid matrix: not a matrix\n");
        return EXIT_FAILURE;
    }
    if (mm_is_dense(mcode) || mm_is_array(mcode) ){
        fprintf(stderr,"invalid matrix: not a sparse matrix\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

/* parse MatrixMarket matrix in @fp into @mat
 * expand simmetric matrix into a normal matrix with both parts
 * return sparse matrix
 */
static int _MMtoCSR(spmat* mat, FILE *fp, MM_typecode mcode) {
    int out=EXIT_FAILURE, scanndRet=0;
    uint nzIdx=0;         //expanded num of nz (in case of sym matrix)
    uint idx;
    double val;         //current entry val
    uint row,col;       //current entry's row,col
    uint* rowsNextCol = NULL;  //for each row -> next entry progressive idx
    int* _rowsLastCol = NULL; //for each row -> last added entry's columnIdx 
    entry* entries = NULL;     //MM parsed and exmpanded entries
    ///consistency checks
    ///init
    if (mm_is_symmetric(mcode)) {
        mat->NZ *= 2;
    }
    if (!(entries         = malloc(mat->NZ * sizeof(*entries))))  goto _free;
    CONSISTENCY_CHECKS{
        if (!(rowsNextCol = calloc(mat->M,sizeof(*rowsNextCol)))) goto _free;
    }
    if (!(_rowsLastCol    = malloc(mat->M*sizeof(*rowsNextCol)))) goto _free;
    memset(_rowsLastCol,-1,mat->M*sizeof(*_rowsLastCol));
    ///parse MM fp lines
    while (1) { // Reading the fp until EOF
        if          (mm_is_pattern(mcode)){
            scanndRet = fscanf(fp, "%u %u\n", &row, &col);
            val = 1.0;
        } else if   (mm_is_real(mcode) || (mm_is_integer(mcode))){
            scanndRet = fscanf(fp, "%u %u %lf\n", &row, &col, &val);
        }
        

        if (scanndRet == EOF){ //TODO more strict check with type&ret?
            if (ferror(fp)){
                perror("fscanf EOF");
                goto _free;
            } else  break;
        }
        CONSISTENCY_CHECKS{ 
            //NNZERO VAL CHECKS 
            if (!val){
                fprintf(stderr,"invalid sparse matrix with 0 entry explicitly stored\n");
                goto _free;
            }
            //TODO USELESS ? ? ?
            if ((mm_is_pattern(mcode) && scanndRet != 2) || 
                (!mm_is_pattern(mcode) && scanndRet != 3)){
                fprintf(stderr,"invalid matrix: not consistent entry scannable\n");
                goto _free;
            }
        }
       ////ADD THE CURRENT MATRIX ENTRY
       mat -> IRP[row]++;   //now just count entries per row
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
    mat->NZ = nzIdx; //all rappresented element (expanding symmetric ones)
    
    if(!(mat->JA = malloc(nzIdx*sizeof(*(mat->JA))))){
        fprintf(stderr,"JA malloc err\n");
        goto _free;
    }
    if(!(mat->AS = malloc(nzIdx*sizeof(*(mat->AS))))){
        fprintf(stderr,"AS malloc err\n");
        goto _free;
    }
    //IRP use rows lens as increments to build row index "pointer"
    //0th -> 0 mandatory; 1th = 0th row len, ...., M+1th = end of Mth row
    for (uint i=2; i<mat->M+1; i++)    mat->IRP[i] += mat->IRP[i-1];
    CONSISTENCY_CHECKS{
        if (mat->IRP[mat->M] != mat->NZ){
            fprintf(stderr,"IRP[M] %u != NZ %u\n",mat->IRP[mat->M],mat->NZ);
            goto _free;
        }
    }
    //build core struct of CSR
    for (uint i=0; i<nzIdx; i++) {
        row = entries[i].row;
        col = entries[i].col;
        val = entries[i].val;
        CONSISTENCY_CHECKS{
            if (_rowsLastCol[row] >= (int) col){
                fprintf(stderr,"not sorted rows\n");
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
        goto out;
    }
    //banner -> parse  matrix specs
    if (mm_read_banner(fp, &mcode) != 0) {
        fprintf(stderr,"mm_read_banner err at:%s\n",matPath);
        goto out;
    }
    //assert matrix is compatible with this app scope
    if (MMCheck(mcode))     goto out;
    //alloc sparse matrix components
    if (!(mat = calloc(1,sizeof(*mat)))){
        fprintf(stderr," mat struct alloc errd");
        goto out;
    }
    //parse sizes
    if(mm_read_mtx_crd_size(fp, &mat->M, &mat->N, &mat->NZ)){
        fprintf(stderr,"mm_read_mtx_crd_size err at %s:\n",matPath);
        goto err;
    }
    if (!(mat->IRP = calloc(mat->M+1,sizeof(*(mat->IRP))))){
        fprintf(stderr,"IRP calloc err\n");
        goto err;
    }

    if (_MMtoCSR(mat, fp, mcode)){
        fprintf(stderr,"MAT PARSE TO CSR ERR at:%s\n",matPath);
        goto err;
    }
    
    out:
    fclose(fp);
    return mat;

    err:
    if(mat->IRP)    free(mat->IRP); 
    if(mat->AS)     free(mat->AS); 
    if(mat->JA)     free(mat->JA); 
    free(mat);
    fclose(fp);
    return NULL;
}
//TODO spmat* MMtoHELL(char* matPath){}
