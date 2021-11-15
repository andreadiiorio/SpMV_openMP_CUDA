#ifndef UTILS
#define UTILS

#include <stddef.h> 
#include "macros.h"

extern int urndFd;	//file pointer to DRNG_DEVFILE O_RDONLY opened
int init_urndfd(); // wrap init urndFd
/*
 * urndFd usage template to populate random timeout
	if(_read_wrap(urndFd,(char*)&timeout,sizeof(timeout))<0){
		fprintf(stderr,"rnd read for thread's timeout failed");
		ret=EXIT_FAILURE;
		goto end;
	}
 */
//wrap read cycle over @fd
int read_wrap(int fd,char* dst,size_t count);

//parse configuration from env
int getConfig(CONFIG* conf);

//append only list implemented with a reallocated array
typedef struct{
    ulong* a;
    ulong  size;
    ulong  lastIdx;
} APPENDARRAY;
//append @val to @list, reallocating if reached end
//TODO inline int appendArr(ulong val,APPENDARRAY* list);

void sortuint(uint* arr, uint len);     //sort uint array @arr of @len elements
void sortulong(ulong* arr, ulong len);   //sort ulong array @arr of @len elements


//return 0 if vectors a and b has elements that differ at most of DOUBLE_DIFF_THREASH 
int doubleVectorsDiff(double* a, double* b, ulong n);
//fill a random vector in @v long @size doubles
int fillRndVector(ulong size, double* v);
//read vector as a sequence of space separated double from file at @fpath 
#define VECTOR_STEP_MALLOC 100

/* 
 * decompress file at @path into @tmpFsDecompressPath, 
 * decompression command obtanined first looking at the extension
 * then matching it with a list of avaible decompression cmd
 * that can be make as shell cmd adding @path > @tmpFsDecompressPath
 * e.g. decompress xz -d -c @path > @tmpFsDecompressPath
 * Returns: -1 if decompression wasn't possible otherwise decompress command exti status
 */
int extractInTmpFS(char* path, char* tmpFsDecompressPath);
//compute E[@values] in @out[0] and VAR[@values] in @out[1] of @numVals values
void statsAvgVar(double* values,uint numVals, double* out);
/*convert @sparseMat sparse matrix in dense matrix returned*/
double* CSRToDense(spmat* sparseMat);
void printMatrix(double* mat,ulong m,ulong n,char justNZMarkers);
void printSparseMatrix(spmat* sparseMat,char justNZMarkers);
void printVector(double* v,ulong size);

#endif
