#ifndef UTILS
#define UTILS

///AUX MACROS
typedef unsigned int uint;
typedef unsigned long ulong;
#define swap(a,b)           a=a^b;b=b^a;a=a^b
#define IDX2D(i,j,nCols)    (j + i*nCols)
#define DRNG_DEVFILE        "/dev/urandom"
///Smart controls
//TODO TOGGLE! ADD
#define TRUE    1
#define FALSE   0
#define DEBUG if( TRUE )
#define CONSISTENCY_CHECKS  if( TRUE )
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
//fill a random vector in @v long @size doubles
int fillRndVector(uint size, double* v);
//read vector as a sequence of space separated double from file at @fpath 
#define VECTOR_STEP_MALLOC 100

/*convert @sparseMat sparse matrix in dense matrix returned*/
double* CSRToDense(spmat* sparseMat);
void printMatrix(double* mat,uint m,uint n,char justNZMarkers);
void printSparseMatrix(spmat* sparseMat,char justNZMarkers);

#endif
