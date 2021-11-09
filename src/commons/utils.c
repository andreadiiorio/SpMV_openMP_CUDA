//various aux functions
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <limits.h>

#include "macros.h"
#include "sparseMatrix.h"
#include "utils.h"


int urndFd; //will point to urandom device file


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

int getConfig(CONFIG* conf){
    int changes=EXIT_FAILURE;
    char *varVal,*ptr;
    ulong val;
    if ((varVal = getenv(GRID_ROWS))){
        val=strtoul(varVal,&ptr,10);
        if (ptr==varVal || val>= UINT_MAX){
            perror("strtol errd");
        } else {
            conf->gridRows = val;
        }
        changes = EXIT_SUCCESS;
    }
    if ((varVal = getenv(GRID_COLS))){
        val=strtoul(varVal,&ptr,10);
        if (ptr==varVal || val>= UINT_MAX){
            perror("strtol errd");
        } else {
            conf->gridCols = val;
        }
        changes = EXIT_SUCCESS;
    }
    return changes;
}
/////LIB-SORTING -- WRAPPERS
//comparing functions
int cmpuint(const void* a, const void*b){
    uint aa=*((uint*) a), bb = *((uint*) b);
    return aa==bb?0:aa>bb?1:-1;
}
//sorting functions 
void sortuint(uint* arr, uint len){
    qsort(arr,len,sizeof(*arr),cmpuint);
}
///MATH UTILS

static inline int rndDouble_sinAll(double* d){
    if(read_wrap(urndFd,(void*) d,sizeof(*d))){
        ERRPRINT("read_wrap failed to read rnd double\n");
        return EXIT_FAILURE;
    }
    *d = sin(*d) * MAXRND;
    return EXIT_SUCCESS;
}
long _rndHold;  //permanent storage of rnd longs
static inline int rndDouble_sinDecimal(double* d){
    if(read_wrap(urndFd,(void*) &_rndHold,sizeof(_rndHold))){
        ERRPRINT("read_wrap failed to read holder for rnd double\n");
        return EXIT_FAILURE;
    }
    *d = (_rndHold % MAXRND) + sin(_rndHold);
    return EXIT_SUCCESS;
}
   
void statsAvgVar(double* values,uint numVals, double* out){
    double sum=0,sumSquare=0;
    for (uint i=0;  i<numVals;  i++){
        sum += values[i];
        sumSquare += values[i]*values[i];
    }
    out[0]  =   sum/numVals;                            //AVG
    out[1]  =   sumSquare/numVals - (out[0] * out[0]);  //VAR
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
    uint nnz = 0;
    for (uint i=0; i<n; i++){
        diff = ABS( a[i] - b[i] );
        if (MAX(a[i],b[i]))     nnz++; //count nnz
        if( diff > DOUBLE_DIFF_THREASH ){
            fprintf(stderr,"DIFF IN DOUBLE VECTORS: %lf > threash=%lf\tat nnz:%u\n",
                diff,DOUBLE_DIFF_THREASH,nnz);
            return EXIT_FAILURE;
        }
        else if (diff > maxDiff)    maxDiff = diff;
    }
    DEBUG{
        printf("checked diff between 2 double vector of %4u nnz with "
          "max diff: %le < threash: %le\n",nnz,maxDiff,DOUBLE_DIFF_THREASH);
        if (!maxDiff){ //self diff check uselss TODO REMOVE
            if (!memcpy(a,b,n*sizeof(*a)))
                printf("exact matching among the 2 double vectors\n!");
        }
    }
    return EXIT_SUCCESS;
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


void printVector(double* v,uint size){
    for( uint i=0;i<size;i++)   printf("%1.1lf ",v[i]);
    printf("\n");
}

////VAR -- MISC

inline int appendArr(uint val,APPENDARRAY* list){
    return 0;   //TODO
}
