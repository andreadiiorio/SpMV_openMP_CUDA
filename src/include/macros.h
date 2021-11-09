#ifndef MACROS
#define MACROS


///aux macro-functions
#define	ABS(a)				        ((a) > 0   ? (a) : -(a))
#define	MIN(a,b)			        ((a) < (b) ? (a) : (b))
#define MAX(a,b)			        ((a) > (b) ? (a) : (b))
#define AVG(a,b)                    ( a/2 + b/2 + (a%2+b%2)/2 )
#define swap(a,b)                   a=a^b;b=b^a;a=a^b
#define MAT_IDX_ROWMAJ(r,c,cols)	( r*cols+c )
//ceil(x/y) with integers
#define INT_DIV_CEIL(x,y)		    ( (x-1) / y + 1 )
//2D ROW MAJOR indexing wrap compute
#define IDX2D(i,j,nCols)            (j + i*nCols)
///distribuite reminder @rem in group givin an extra +1 to the first @rem
#define UNIF_REMINDER_DISTRI(i,div,rem) \
    ( div + ( i<rem ? 1 : 0 ) )
#define UNIF_REMINDER_DISTRI_STARTIDX(i,div,rem) \
    ( i * div + MIN(i,rem)*1 )

////PRINTS
#define CHIGHLIGHT                  "\33[1m\33[92m"
#define CHIGHLIGHTERR               "\33[31m\33[1m\33[44m"
#define CEND                        "\33[0m"
#define hprintsf(str,...)           printf( CHIGHLIGHT str CEND,__VA_ARGS__ ) 
#define hprintf(str)                printf( CHIGHLIGHT str CEND) 
#define ERRPRINTS(str,...)          fprintf( stderr, CHIGHLIGHTERR str CEND,__VA_ARGS__ )
#define ERRPRINT(str)               fprintf( stderr, CHIGHLIGHTERR str CEND )
///CONSTANTS
#define DOUBLE_DIFF_THREASH         7e-5
#define DRNG_DEVFILE                "/dev/urandom"
#define MAXRND                      1996

///Smart controls
#define FALSE                       ( 0 )
#define TRUE                        ( ! FALSE )
//debug checks and tmp stuff
#ifndef DEBUG 
    #define DEBUG                       if( TRUE )
#endif
//long prints
#ifndef DEBUGPRINT
    #define DEBUGPRINT                  if( FALSE )
#endif
//heavy impact debug checks
#ifndef DEBUGCHECKS
    #define DEBUGCHECKS                 if( FALSE )
#endif
//extra print in the normal output
#ifndef AUDIT_INTERNAL_TIMES
    #define AUDIT_INTERNAL_TIMES        if( TRUE )
#endif
#ifndef VERBOSE
    #define VERBOSE                     if( TRUE )
#endif
//extra checks over the imput and correct partials
#ifndef CONSISTENCY_CHECKS
    #define CONSISTENCY_CHECKS          if( TRUE )
#endif
#ifndef AVG_TIMES_ITERATION
    #define AVG_TIMES_ITERATION         5
#endif
///aux types
typedef unsigned char  uchar;
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


//////TODO CONFIGURATION DEFINITIONS
extern double Start,End,Elapsed,ElapsedInternal;
typedef struct{
    ushort gridRows;
    ushort gridCols;
    //TODO FULL CONFIG DOCCED HERE
    uint threadNum;  //thread num to use in an OMP parallel region ...
    void* spgemmFunc;   //aux spgemm function to use. 
    //TODO MAKE THIS A CONTAINER OF SUB STRUCT PASSABLE TO SPGEMM FUNCS TO AVOID CAST
} CONFIG;  
///config from ENV
#define GRID_ROWS   "GRID_ROWS"
#define GRID_COLS   "GRID_COLS"
#endif 	//MACROS
