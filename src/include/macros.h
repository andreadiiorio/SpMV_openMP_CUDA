#ifndef MACROS
#define MACROS


///aux macro-functions
#define	ABS(a)				((a) > 0   ? (a) : -(a))
#define	MIN(a,b)			((a) < (b) ? (a) : (b))
//#define MAX(a,b)			((a) > (b) ? (a) : (b))
#define swap(a,b)           a=a^b;b=b^a;a=a^b
#define MAT_IDX_ROWMAJ(r,c,cols)	( r*cols+c )
//ceil(x/y)
#define INT_CEIL_DIV(x,y)		( (x-1) / y + 1 )
#define IDX2D(i,j,nCols)    (j + i*nCols)

#define ERRPRINT(str)       fprintf(stderr,str)
///CONSTANTS
#define DOUBLE_DIFF_THREASH 1e-3
#define DRNG_DEVFILE        "/dev/urandom"
#define MAXRND              1996
///Smart controls
//TODO TOGGLE! ADD
#define FALSE               ( 0 )
#define TRUE                ( ! FALSE )
#define DEBUG               if( TRUE )
#define VERBOSE             if( TRUE )
#define CONSISTENCY_CHECKS  if( TRUE )
///aux types
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long  ulong;
//smart decimal type custom precision def build macro _DECIMAL_TRGT_PREC 
#ifndef _DECIMAL_TRGT_PREC
//dflt floating point precision & formatting chars
	#define _DECIMAL_TRGT_PREC	double
	#define _DECIMAL_TRGT_PREC_PR 	"%lf"
#else 
    //TODO SELECT WITH RESPECT TO THE EXPORTED TARGET DECIMAL TYPE
	#define _DECIMAL_TRGT_PREC_PR 	"%f"
#endif
typedef _DECIMAL_TRGT_PREC	decimal;


///EXTRA INCLUDE    --- cuda 
///assertionn are disabled at compile time by defining the NDEBUG preprocessor macro before including assert.h	s
#ifdef ASSERT
	#include <assert.h>
#endif


#endif 	//MACROS
