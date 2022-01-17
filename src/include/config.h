#ifndef CONFIG_H
#define CONFIG_H

typedef struct{
    ushort gridRows;
    ushort gridCols;
    //TODO FULL CONFIG DOCCED HERE
    uint   threadNum;  //thread num to use in an OMP parallel region ...
    void*  chunkDistrbFunc;  //CHUNKS_DISTR_INTERF func pntr
	#ifdef __CUDACC__
	dim3   gridSize;
	dim3   blockSize;
	size_t sharedMemSize;
	#endif
} CONFIG;  
///Smart controls
#define FALSE                       ( 0 )
#define TRUE                        ( ! FALSE )
///AUDIT&CHECKS
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
    #define VERBOSE                     if( FALSE )
#endif
//extra checks over the imput and correct partials
#ifndef CONSISTENCY_CHECKS
    #define CONSISTENCY_CHECKS          if( TRUE )
#endif
///AUDIT extra configuration
//#define ROWLENS
#ifdef ROWLENS
#pragma message("ROW_LENS ARRAY ENABLED")
#endif
///CONSTS
#define ELL_MAX_ENTRIES ( 100l << 20 )  //ell max (padded) entries supported in a matrix 
#define ELL_AS_FILLER       (0 )        //handled with calloc
//TODO NOW FILLED WITH LAST NNPADDED COL #define ELL_JA_FILLER       (-1)
//#define DOUBLE_VECT_DIFF_EARLY_EXIT 1
//#define RNDVECTMIN          222222
#define VECTOR_STEP_REALLOC 25
#define VECTOR_READ_BLOCK	50		//file (raw) vector read block
#define RNDVECTORSIZE       100000
#define RNDVECTORDUMP       TMPDIR  "rndVectorDump"
#define RNDVECTORDUMPRAW    TMPDIR  "rndVectorDumpRaw"
#define OUTVECTORDUMP       TMPDIR  "outVectorDump"
#define OUTVECTORDUMPRAW    TMPDIR  "outVectorDumpRaw"
//#define FLOAT_DIFF_ABS
#ifndef AVG_TIMES_ITERATION
    #define AVG_TIMES_ITERATION         5
#endif
//ompChunksDivide.h -> chunksFairFolded()
#ifndef FAIR_CHUNKS_FOLDING
    #define FAIR_CHUNKS_FOLDING 4
#endif
//SPMV specific
//rows partitions for dotProduct SIMD reduction enable
#ifndef SIMD_ROWS_REDUCTION
    #define SIMD_ROWS_REDUCTION         TRUE
#endif
#if SIMD_ROWS_REDUCTION == TRUE
    #pragma message("SIMD_ROWS_REDUCTION enabled")
    //TODO SOME TRICK TO HAVE 1! PRINT
#endif
extern double Start,End,Elapsed,ElapsedInternal;
#define DOUBLE_DIFF_THREASH         7e-5
#define DRNG_DEVFILE                "/dev/urandom"
#define MAXRND                      1
#ifndef TMPDIR
    #define TMPDIR                      "/tmp/"
#endif
#define TMP_EXTRACTED_MARTIX    TMPDIR "extractedMatrix"

#endif  //CONFIG_H
