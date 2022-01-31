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

#ifndef UTILS
#define UTILS

#include <stddef.h> 
#include "macros.h"

extern int urndFd;	//file pointer to DRNG_DEVFILE O_RDONLY opened
int init_urndfd(); // wrap init urndFd

///IO
//UNBUFFERED IO
/*
 * urndFd usage template to populate random timeout
	if(_read_wrap(urndFd,(char*)&timeout,sizeof(timeout))<0){
		fprintf(stderr,"rnd read for thread's timeout failed");
		ret=EXIT_FAILURE;
		goto end;
	}
 */
//wrap read cycle over @fd
int read_wrap(int fd,void* dst,size_t count);
int fread_wrap(FILE* fp,void* dst,size_t count);
//dual of read_wrap
int write_wrap(int fd,void* src,size_t count);
//create or open file at @outFpath for write
int createNewFile(char* const outFpath);
///STRUCTURED DATA IO
#define DOUBLE_STR_FORMAT	"%25le\n"
//write double vector @v as row sequence of double at @fpath
//e.g. read with od -tf8 -w8 fpath : OCTALOFFSET:   DOUBLE FULL DIGITS
int writeDoubleVector(char* fpath,double* v,ulong size);
/*
 * read vector of double [Str] of arbitrary size from @fpath, true lenght in *size
 * if size point to a nnz value, the initial allocation will be of *size
 * eventual successive reallocation done multipling *size with VECTOR_STEP_REALLOC
 */
double* readDoubleVector(char* fpath,ulong* size);
double* readDoubleVectorStr(char* fpath,ulong* size);

///STRUCTURED DATA IO -- BUFFERED: FSCANF - FPRINTF
//dual of readDoubleVectorVector
int writeDoubleVectorAsStr(char* fpath,double* v,ulong size);

#include "config.h"
///config from ENV
#define GRID_ROWS   "GRID_ROWS"
#define GRID_COLS   "GRID_COLS"
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


/*
 * return 0 if vectors a and b has elements that differ at most of DOUBLE_DIFF_THREASH 
 * if diffMax!=NULL save there the max difference value  
 *  of the 2 entries of the input vectors, signed as a[i] - b[i] (as well as dump prints)
 * CONVENTION:	@a = true result, @b vector to check with the true result 
 */
int doubleVectorsDiff(double* a, double* b, ulong n,double* diffMax);
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
