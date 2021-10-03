//various aux functions
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "sparseMatrix.h"
#include "utils.h"

int urndFd; //from utils

//rnd gen from /dev/random
int read_wrap(int fd,char* dst,size_t count){
	ssize_t rd;
	size_t readed=0;
	while (readed < count){
		rd=read(fd,dst+readed,count-readed);
		if (rd<0){
			perror("read");
			return rd;
		}
		readed+=rd;
	}
	return 0;
}

int init_urndfd(){ // wrap init urndFd
	if((urndFd=open(DRNG_DEVFILE,O_RDONLY))<0){
			perror("open DRNG_DEVFILE");
			return EXIT_FAILURE;
	}
    return EXIT_SUCCESS;
}

/// MATRIX - VECTOR UTILS
int fillRndVector(uint size, double* v){
    for( uint x=0; x<size; ++x ){
        if(read_wrap(urndFd,(char*)v+x,sizeof(*v))<0){
	    	fprintf(stderr,"rnd read for thread's timeout failed");
	    	return EXIT_FAILURE;
	    }
    }
    return EXIT_SUCCESS;
}

double* CSRToDense(spmat* sparseMat){
    double* denseMat;
    uint i,j,idxSparse;
    if (!(denseMat = calloc(sparseMat->M*sparseMat->N, sizeof(*denseMat)))){
        fprintf(stderr,"dense matrix alloc failed\n");
        return  NULL;
    }
    for (i=0;i<sparseMat->M;i++){
        for (idxSparse=sparseMat->IRP[i];idxSparse<sparseMat->IRP[i+1];++idxSparse){
             j = sparseMat->JA[idxSparse];
             //converting sparse item into dense entry
             denseMat[IDX2D(i,j,sparseMat->N)] = sparseMat->AS[idxSparse]; 
        }
    }
    return denseMat;
}

void printMatrix(double* mat,uint m,uint n,char justNZMarkers){
    printf("printing matrix: %u x %u\n",m,n);
    uint i,j;
    for (i=0;i<m;i++){
        for (j=0;j<n;j++){
            if (justNZMarkers)  printf("%s",mat[IDX2D(i,j,n)]?".":" ");
            else                printf("%1.1lf ",mat[IDX2D(i,j,n)]);
        }
        printf("\n");
    }
}

void printSparseMatrix(spmat* spMatrix,char justNZMarkers){
    double* denseMat = CSRToDense(spMatrix);
    if (!denseMat)  return;
    printMatrix(denseMat,spMatrix->M,spMatrix->N,justNZMarkers);
    free(denseMat);
    }
