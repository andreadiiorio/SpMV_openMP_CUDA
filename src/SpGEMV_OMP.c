#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "SpGEMV.h"
#include "parser.h"
#include "utils.h"
#include "sparseMatrix.h"

int main(int argc, char** argv){
    int out=EXIT_FAILURE;
    if (init_urndfd())  return out;
    if (argc < 3 ){
        fprintf(stderr,
        "usage: MatrixMarket sparse matrix COO, vectorFile || "RNDVECT"\n");
        return out;
    }
    double *vector = NULL, *outVector = NULL;
    spmat* mat = NULL; 
    ///parse sparse matrix and dense vector
    if (!(mat = MMtoCSR(argv[1]))){
        fprintf(stderr,"err during conversion MM -> CSR\n");
        return out;
    }
    //get the vector
    if (!(strncmp(argv[2],RNDVECT,strlen(RNDVECT)))){
        //generate a random vector
        if (!(vector = malloc(RNDVECTORSIZE * sizeof(*vector)))){
            fprintf(stderr,"rnd vector malloc failed\n");
            goto _free;
        }
        if (fillRndVector(RNDVECTORSIZE,vector)){
            fprintf(stderr,"fillRndVector errd\n");
            goto _free;
        }
    } else{
        //read vector from the given file
        if (!(vector = readVector(argv[2]))){
            fprintf(stderr,"err during readVector at:%s\n",argv[2]);
            goto _free;
        }
    }
    if (!(outVector = malloc( mat->M * sizeof(*outVector)))){
        fprintf(stderr,"outVector malloc errd\n");
        goto _free;
    }
    
    DEBUG        printSparseMatrix(mat,TRUE);
    ///TODO PARALLEL COMPUTATION
        

    out = EXIT_SUCCESS;
    _free:
    if (mat)          free(mat);
    if (vector)       free(vector);
    if (outVector)    free(outVector);
    return out;
}
