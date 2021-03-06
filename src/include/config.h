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
	#pragma message("DEBUG MACROS ENABLED")
    #define DEBUG                       if( TRUE )
#endif
//long prints
#ifndef DEBUGPRINT
	#pragma message("DEBUGPRINT MACROS ENABLED !!!")
    #define DEBUGPRINT                  if( FALSE )
#endif
//heavy impact debug checks
#ifndef DEBUGCHECKS
	#pragma message("DEBUGCHECKS MACROS ENABLED !!!")
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
#define ELL_MAX_ENTRIES ( 6l << 27 ) //2*6GB of ell (padded) entries maxSupport in a matrix 
#define LIMIT_ELL_SIZE				 //enable above threshold
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

//#ifdef __CUDACC_	//CUDA CONFIG //TODO not defined here
///threading sizeing
#ifndef BLOCKS_1D
	#define BLOCKS_1D			( 1u << 8 )
	#pragma message("CUDA BLOCKS_1D: " STRIFY(BLOCKS_1D))
#endif
#ifndef BLOCKS_2D_WARP_R
	#define BLOCKS_2D_WARP_R	( 1u << (10-5) )
	#pragma message("CUDA BLOCKS_2D_WARP_R: " STRIFY(BLOCKS_2D_WARP_R))
#endif
//#endif	//__CUDACC_

extern double Start,End,Elapsed,ElapsedInternal;
#define DOUBLE_DIFF_THREASH         7e-4
#define DRNG_DEVFILE                "/dev/urandom"
#define MAXRND                      3e-5
#ifndef TMPDIR
    #define TMPDIR                      "/tmp/"
#endif
#define TMP_EXTRACTED_MARTIX    TMPDIR "extractedMatrix"

#endif  //CONFIG_H
