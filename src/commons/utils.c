//various aux functions
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>

#include "sparseMatrix.h"
#include "utils.h"
#include "macros.h"

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

int readALL_wrap(int fd,char** dst,size_t* count){
	ssize_t rd=!0; //to allow count > fsize
	size_t readed=0;
    char allocated=0;   //flag if required *dst allocation
    if (!(*dst)){  //allocate dst buffer of same size of file
        off_t seekCurr=lseek(fd,0,SEEK_CUR);
        off_t fsize=lseek(fd,0,SEEK_END);
        if( seekCurr==-1 || fsize==-1 || lseek(fd,seekCurr,SEEK_SET)==-1){
            perror("lseek");
            return EXIT_FAILURE;
        }
        *count=fsize;
        if (!(*dst=malloc(fsize))){
            fprintf(stderr,"malloc read_wrap file size buf error\n");
            return EXIT_FAILURE;
        }
        allocated=!0;
    }
    //read loop
	while (readed < *count && rd > 0){
		rd=read(fd,(*dst)+readed,*count-readed);
		if (rd<0){
			perror("read");
            if (allocated) free(*dst);
			return rd;
		}
		readed+=rd;
	}
    if (readed < *count) (*dst)[readed]='\0';    //TODO NEEDED?
	return EXIT_SUCCESS;
}
int init_urndfd(){ // wrap init urndFd
	if((urndFd=open(DRNG_DEVFILE,O_RDONLY))<0){
			perror("open DRNG_DEVFILE");
			return EXIT_FAILURE;
	}
    return EXIT_SUCCESS;
}
int createNewFile(char* const outFpath){
    int mode=S_IRWXU;
    int outFd=open(outFpath, O_WRONLY | O_CREAT | O_TRUNC, mode);
    //TODO ? RIMETTI O_EXCL if (errno==EEXIST)      outFd=open(outFpath, O_WRONLY | O_TRUNC, mode);
    if (outFd<0)            perror("open outFd failed ");
    return outFd;
}

///MATH UTILS

inline int rndDouble_sinAll(double* d){
    if(read_wrap(urndFd,(void*) d,sizeof(*d))){
        ERRPRINT("read_wrap failed to read rnd double\n");
        return EXIT_FAILURE;
    }
    *d = sin(*d) * MAXRND;
    return EXIT_SUCCESS;
}
long _rndHold;  //permanent storage of rnd longs
inline int rndDouble_sinDecimal(double* d){
    if(read_wrap(urndFd,(void*) &_rndHold,sizeof(_rndHold))){
        ERRPRINT("read_wrap failed to read holder for rnd double\n");
        return EXIT_FAILURE;
    }
    *d = (_rndHold % MAXRND) + sin(_rndHold);
    return EXIT_SUCCESS;
}
     

/// MATRIX - VECTOR UTILS
int fillRndVector(uint size, double* v){
    for( uint x=0; x<size; ++x ){
        if(rndDouble_sinDecimal( v+x )) return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int doubleVectorsDiff(double* a, double* b, uint n){
    double diff,maxDiff=0;
    for (uint i=0; i<n; i++){
        diff = ABS( a[i] - b[i] );
        if( diff > DOUBLE_DIFF_THREASH ){
            fprintf(stderr,"DIFF IN DOUBLE VECTORS: %lf > threash: %lf",
                diff,DOUBLE_DIFF_THREASH);
            return EXIT_FAILURE;
        }
        else if (diff > maxDiff)    maxDiff = diff;
    }
    VERBOSE printf("checked diff between 2 double vector with "
        "max diff: %lf < threash: %lf\n",maxDiff,DOUBLE_DIFF_THREASH);
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
    printf("\n");
}

void printSparseMatrix(spmat* spMatrix,char justNZMarkers){
    double* denseMat = CSRToDense(spMatrix);
    if (!denseMat)  return;
#ifdef ROWLENS
    for (ushort i=0; i<spMatrix->M; i++)    printf("%u\t%u\n",i,spMatrix->RL[i]);
#endif
    printMatrix(denseMat,spMatrix->M,spMatrix->N,justNZMarkers);
    free(denseMat);
}

void printVector(double* v,uint size){
    for( uint i=0;i<size;i++)   printf("%1.1lf ",v[i]);
    printf("\n");
}
