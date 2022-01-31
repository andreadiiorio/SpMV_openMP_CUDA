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

#ifndef MACROS
#define MACROS


///aux macro-functions
#define	ABS(a)				        ((a) > 0   ? (a) : -(a))
#define	MIN(a,b)			        ((a) < (b) ? (a) : (b))
#define MAX(a,b)			        ((a) > (b) ? (a) : (b))
#define AVG(a,b)                    ((a)/2 + (b)/2 + ((a)%2+(b)%2)/2)
#define SWAP(a,b)                   (a)=(a)^(b);(b)=(b)^(a);(a)=(a)^(b)
//ceil(x/y) with integers
#define INT_DIV_CEIL(x,y)		    ( ( (x) - 1) / (y) + 1 )
//2D ROW MAJOR indexing wrap compute
#define IDX2D(i,j,nCols)            ((j) + (i)*(nCols))
///distribuite reminder @rem in group givin an extra +1 to the first @rem
#define UNIF_REMINDER_DISTRI(i,div,rem) \
    ( (div) + ( (i) < (rem) ? 1 : 0 ) )
#define UNIF_REMINDER_DISTRI_STARTIDX(i,div,rem) \
    ( (i) * (div) + MIN( (i),(rem) ) )

#define STATIC_ARR_ELEMENTS_N(arr)  (sizeof( (arr) ) / (sizeof(*(arr))))  
///STR UTILS
#define _STRIFY(x)  #x
#define STRIFY(x)   _STRIFY(x)
//TODO ADD IN THE COMPARE THE LAST \0 FOR FULL STRING MATCH AND NOT PARTIAL
#define strEqual(s0,s1)	!(strncmp( (s0) , (s1) , strlen( (s1) )))	
	//convention of putting the target->template as s1

#include <stdlib.h>
#include <stdio.h>

////PRINTS
#define CHIGHLIGHT                  "\33[1m\33[92m"
#define CCC                         CHIGHLIGHT
#define CHIGHLIGHTERR               "\33[31m\33[1m\33[44m"
#define CCCERR                      CHIGHLIGHTERR
#define CEND                        "\33[0m"
#define hprintsf(str,...)           printf( CHIGHLIGHT str CEND,__VA_ARGS__ ) 
#define hprintf(str)                printf( CHIGHLIGHT str CEND) 
#define ERRPRINTS(str,...)          fprintf( stderr, CHIGHLIGHTERR str CEND,__VA_ARGS__ )
#define ERRPRINT(str)               fprintf( stderr, CHIGHLIGHTERR str CEND )

#include <assert.h> 

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
//#ifdef ASSERT 	#include <assert.h> #endif

#endif 	//MACROS
