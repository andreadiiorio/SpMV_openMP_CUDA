Matrix Market: C routines for Matrix Market I/O
<!--
function enist(u){
 var q="'",c=":",d=".",a="&#64;",s="/",e="=",b=" ";
 i="ilto",t="a",m="ma",h="href",g="gov",n="nist";
 var x=u+a+n+d+g,l="<",g=">";
 document.write(l+t+b+h+e+q+m+i+c+x+q+g+x+l+s+t+g); }
//-->
### Matrix Market
ANSI C library for Matrix Market I/O
====================================
The numerical data in the Matrix Market file formats can be easily processed
using variants of **fscanf()** and **fprintf()**
functions.
The only non-trivial issue is to figure out what kind of matrix is
represented in a Matrix Market file. Because of the wide range of
possibilities, it is impractical to have a single function handle
every case (furthermore, most applications will support only a subset
of these matrix types). Instead

__we provide utilities that identify and manage only the type and size information in MM files,
leaving the actual reading and writing mechanisms to the driving application or higher-level I/O routines
   -> one matrix entry per line:    <ROWID COLID VALUE>__

__Reading a Matrix Market into three basic steps__
1. use [mm\_read\_banner()](#mm_read_banner) to process the 1st line of file and identify the matrix type
2. use a type-specific function, such as
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size)
 to skip the optional comments and process the matrix size information
3. use a variant of scanf() to read the numerical data, one matrix entry per line

__ Saving a matrix from an internal structure to Matrix Market __
1. use [mm\_write\_banner()](#mm_write_banner)
to create the 1st line of
 the Matrix Market file
2. (optional) add '%' delimited comments
3. use a type-specific function, such as
[mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size)
to record the matrix size information
4. use a variant of  printf() to write the numerical data, one matrix entry per line

=========================== CORE API =============================================
##CORE API
mm_read_banner (FILE\*,MM_typecode) -> set  \*MM_typecode<-MMfileInfos
mm_write_banner(FILE\*,MM_typecode)  -> read MM_typecode -> write to MM file header 
##TYPECODE -> HEADER INFO STRUCT
typedef char MM_typecode[4];

##internal definitions: __typecode__

   MM_matrix_typecode: 4-character sequence

				    ojbect 		sparse/   	data        storage 
						  		dense     	type        scheme

   string position:	 [0]        [1]			[2]         [3]

   Matrix typecode:  M(atrix)  C(oord)		R(eal)   	G(eneral)
						        A(array)	C(omplex)   H(ermitian)
											P(attern)   S(ymmetric)
								    		I(nteger)	K(kew)

##QUERY (SET==> is -->set)  OSS typecode "var" in macros triv.ref to MM_typecode
#define mm_is_matrix(typecode)	    ((typecode)[0]=='M')

#define mm_is_sparse(typecode)	    ((typecode)[1]=='C')
#define mm_is_coordinate(typecode)  ((typecode)[1]=='C')
#define mm_is_dense(typecode)	    ((typecode)[1]=='A')
#define mm_is_array(typecode)	    ((typecode)[1]=='A')

#define mm_is_complex(typecode)	    ((typecode)[2]=='C')
#define mm_is_real(typecode)		((typecode)[2]=='R')
#define mm_is_pattern(typecode)	    ((typecode)[2]=='P')
#define mm_is_integer(typecode)     ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)   ((typecode)[3]=='S')
#define mm_is_general(typecode)	    ((typecode)[3]=='G')
#define mm_is_skew(typecode)	    ((typecode)[3]=='K')
#define mm_is_hermitian(typecode)   ((typecode)[3]=='H')

========================================================================
###  Source code
* examples:
	+ [example\_read.c](http://math.nist.gov/MatrixMarket/mmio/c/example_read.c)+ [example\_write.c](http://math.nist.gov/MatrixMarket/mmio/c/example_write.c)* library routines:
	+ [mmio.h](http://math.nist.gov/MatrixMarket/mmio/c/mmio.h)+ [mmio.c](http://math.nist.gov/MatrixMarket/mmio/c/mmio.c)
---
### Documentation
* [mm\_read\_banner()](#mm_read_banner)* [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size)* [mm\_read\_mtx\_array\_size()](#mm_read_mtx_array_size)* [mm\_write\_banner()](#mm_write_banner)* [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size)* [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size)
---
mm\_read\_banner()
------------------
### NAME
mm\_read\_banner - determine the type of matrix being represented in a
 Matrix Market file
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	int mm\_read\_banner(FILE *f, MM\_typecode *t);
```
### DESCRIPTION
mm\_read\_banner() processes the the first line of a Matrix Market file,
e.g.   
 %%MatrixMarket matrix coordinate real general
  
and returns the matrix characteristics.
File descriptor ***f*** is defined in "stdio.h" and is assumed
to have been opened for read access. The predefined descriptor
***stdin*** can be used to read from the standard
input.
***t*** points to an internal structure that describes the
 matrix charateristics. This MM\_typecode is
 more efficient and convenient
 than storing the explicit banner.
 [query functions](#query), such as
 mm\_is\_complex(t), are available to extract this information.
### RETURN VALUES
mm\_read\_banner() returns 0 if succesful. Otherwise,
 it returns one of the error codes below.
###  ERRORS
* MM\_PREMATURE\_EOF if all items are not present on first line of file.
 * MM\_NO\_HEADER if the file does not begin with "%%MatrixMarket".
 * MM\_UNSUPPORTED\_TYPE if not recongizable description.
###  EXAMPLES
 See [example\_read.c](http://math.nist.gov/MatrixMarket/mmio/c/example_read.c).
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_array_size),
 [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size)
  
  
---
mm\_read\_mtx\_crd\_size()
---------------------------
### NAME
mm\_read\_mtx\_crd\_size - read the size information of a sparse matrix
 (coordinate format) in a Matrix Market file
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	int mm\_read\_mtx\_crd\_size(FILE *f, int *M, int *N, int *nz);
```
### DESCRIPTION
After processing the Matrix Market [banner](#mm_read_banner),
mm\_read\_mtx\_crd\_size() reads
past the optional comments and initalizes size variables ***M***
(number of rows), ***N*** (number of columns), and
***nz***(number of non-zeros).
It is assumed that the matrix being read is in coordinate format.
###  RETURN VALUES
mm\_read\_mtx\_crd\_size() return 0 is successful, otherwise
 it returns one of the error codes below.
###  ERRORS
* MM\_PREMATURE\_EOF if an end-of-file is encountered before processing
 these three values.
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_crd_size),
 [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size),
  
  
---
mm\_read\_mtx\_array\_size()
-----------------------------
### NAME
mm\_read\_mtx\_array\_size - read the size information of a dense matrix
 (array format) in a Matrix Market file
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	int mm\_read\_mtx\_array\_size(FILE *f, int *M, int *N);
```
### DESCRIPTION
After processing the Matrix Market [banner](#mm_read_banner),
mm\_read\_mtx\_crd\_size() reads
past the optional comments and initalizes matrix size variables
***M*** (number of rows), ***N*** (number of columns).
It is assumed that the matrix being read is in array format.
###  RETURN VALUES
mm\_read\_mtx\_array\_size() return 0 is successful, otherwise
 it returns one of the error codes below.
###  ERRORS
* MM\_PREMATURE\_EOF if an end-of-file is encountered before reading
 M and N.
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_crd_size),
 [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size),
  
  
---
mm\_write\_banner()
-------------------
### NAME
mm\_write\_banner - record matrix type information in Matrix Market file
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	int mm\_write\_banner(FILE *f, MM\_typecode *t);
```
### DESCRIPTION
mm\_write\_banner() prints the first line of a Matrix Market file,
which consists of the "%%MatrixMarket" token followed by an attribute
list,
e.g.   
 %%MatrixMarket matrix coordinate real general
  
File descriptor ***f*** is defined in "stdio.h" and is assumed
to have been opened for write access. The predefined descriptor
***stdout*** can be used to read from the standard
output.
***t*** points to an internal structure that describes the
 matrix charateristics. This MM\_typecode is
 more efficient and convenient
 than storing the explicit banner. Various
 [assign functions](#query), such as
 mm\_set\_complex(&t), can be used to set
 these characterisitcs.
### RETURN VALUES
mm\_write\_banner() returns 0 if succesful. Otherwise,
 it returns MM\_COULD\_NOT\_OPEN\_WRITE\_FILE.
###  EXAMPLES
 See [example\_write.c](http://math.nist.gov/MatrixMarket/mmio/c/example_write.c).
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_array_size),
 [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size)
  
  
---
mm\_write\_mtx\_crd\_size()
----------------------------
### NAME
mm\_write\_mtx\_crd\_size - write the size information of a dense matrix
 (array format) to a Matrix Market file
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	int mm\_write\_mtx\_crd\_size(FILE *f, int M, int N, int nz);
```
### DESCRIPTION
Record the matrix dimensions ***M*** x ***N*** and
total number of nonzeros, ***nz*** to Matrix Market file.
Typically called after
[mm\_write\_banner()](#mm_write_banner).
###  RETURN VALUES
mm\_write\_mtx\_crd\_size() returns 0 is successful, otherwise
 MM\_COULD\_NOT\_WRITE\_FILE.
###  DIAGNOSTICS
This is a trivial function to write three integers to ***f***
using fprintf().
It is included in the library only for the sake of completeness,
as a counterpart to [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size).
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_crd_size),
 [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size),
  
  
---
mm\_write\_mtx\_array\_size()
------------------------------
### NAME
mm\_write\_mtx\_array\_size - write the size information of a dense matrix
 (array format) to a Matrix Market file
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	int mm\_write\_mtx\_array\_size(FILE *f, int M, int N);
```
### DESCRIPTION
Record the matrix dimensions ***M*** x ***N***
to Matrix Market file.
Typically called after
[mm\_write\_banner()](#mm_write_banner).
###  RETURN VALUES
mm\_write\_mtx\_array\_size() returns 0 is successful, otherwise
 MM\_COULD\_NOT\_WRITE\_FILE.
### DIAGNOSTICS
This is a trivial function to write two integers to ***f***
using fprintf().
It is included in the library only as a counterpart to
[mm\_read\_mtx\_array\_size()](#mm_read_mtx_array_size).
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_crd_size),
 [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size),
  
  
---
MM\_IS
-------
### NAME
mm\_is\_matrix,
mm\_is\_sparse,
mm\_is\_coordinate,
mm\_is\_dense,
mm\_is\_array,
mm\_is\_complex,
mm\_is\_real ,
mm\_is\_pattern,
mm\_is\_integer,
mm\_is\_symmetric,
mm\_is\_general,
mm\_is\_skew,
mm\_is\_hermitian -matrix type query functions
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	int mm\_is\_matrix(MM\_typecode t);
	int mm\_is\_sparse(MM\_typecode t);
	int mm\_is\_coordinate(MM\_typecode t);
	int mm\_is\_dense(MM\_typecode t);
	int mm\_is\_array(MM\_typecode t);
	int mm\_is\_complex(MM\_typecode t);
	int mm\_is\_real(MM\_typecode t);
	int mm\_is\_pattern(MM\_typecode t);
	int mm\_is\_integer(MM\_typecode t);
	int mm\_is\_symmetric(MM\_typecode t);
	int mm\_is\_general(MM\_typecode t);
	int mm\_is\_skew(MM\_typecode t);
	int mm\_is\_hermitian(MM\_typecode t);
```
### DESCRIPTION
MM\_QUERY functions provide a boolean test on matrix typecodes
(MM\_typecode) variables for a given storage property.
The functions return 0 if false, 1 otherwise. Note that these properties
refer only to the storage scheme, not the mathematical properties.
For example, a mathematically symmetric matrix stored with both upper
and lower triangular halves will evaluate mm\_is\_symmetric() false.
###  RETURN VALUES
 MM\_QUERY functions return 1 if true, 0 if false.
###  DIAGNOSTICS
These functions are implemented as macros.
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_array_size),
 [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size),
  
  
---
MM\_SET
--------
### NAME
mm\_set\_matrix,
mm\_set\_sparse,
mm\_set\_coordinate,
mm\_set\_dense,
mm\_set\_array,
mm\_set\_complex,
mm\_set\_real ,
mm\_set\_pattern,
mm\_set\_integer,
mm\_set\_symmetric,
mm\_set\_general,
mm\_set\_skew,
mm\_set\_hermitian,
mm\_clear\_typecode,
mm\_initialize\_typecode --matrix type query functions
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	int mm\_set\_matrix(MM\_typecode &t);
	int mm\_set\_sparse(MM\_typecode &t);
	int mm\_set\_coordinate(MM\_typecode &t);
	int mm\_set\_dense(MM\_typecode &t);
	int mm\_set\_array(MM\_typecode &t);
	int mm\_set\_complex(MM\_typecode &t);
	int mm\_set\_real(MM\_typecode &t);
	int mm\_set\_pattern(MM\_typecode &t);
	int mm\_set\_integer(MM\_typecode &t);
	int mm\_set\_symmetric(MM\_typecode &t);
	int mm\_set\_general(MM\_typecode &t);
	int mm\_set\_skew(MM\_typecode &t);
	int mm\_set\_hermitian(MM\_typecode &t);
	int mm\_clear\_typecode(MM\_typecode &t);
	int mm\_initialize\_typecode(MM\_typecode &t);
```
### DESCRIPTION
MM\_SET functions are used to encode the matrix charateristics to
a Matrix Market file, in conjunction with [mm\_write\_banner()](#mm_write_banner). To use properly, ***t*** first
must be initialized. For exmaple, the set a MM\_typecode to
describe a complex, Hermititan, sparse matrix, one can write
```
	mm\_initialize\_typecode(&t);
	mm\_set\_matrix(&t);
	mm\_set\_complex(&t);
	mm\_set\_coordinate(&t);
	mm\_set\_hermitian(&t);
```
***t*** can then be used as an argument to
mm\_write\_banner().
###  RETURN VALUES
 MM\_SET functions return void.
###  DIAGNOSTICS
These functions are implemented as macros.
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_array_size),
 [mm\_write\_mtx\_crd\_size()](#mm_write_mtx_crd_size),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size),
  
  
---
mm\_typecode\_to\_str()
------------------------
### NAME
mm\_typecode\_to\_str - convert a typecode to a descriptive string
### SYNPOSIS
```
	#include <stdio.h>
	#include "mmio.h"
	char *mm\_typecode\_to\_str(MM\_typecode t);
```
### DESCRIPTION
mm\_typecode\_to\_str converts a MM\_typecode 
to mnemonic string. This can be used as part of a diagnostics
and error reporting. (See "[example\_read.c](http://math.nist.gov/MatrixMarket/mmio/c/example_read.c)"
for possible use.)
###  RETURN VALUES
mm\_typecode\_to\_str() returns a new character string if
 successful, or NULL if an error occured.
###  DIAGNOSTICS
mm\_typecode\_to\_str uses strdup() internally, and
consequently the returned string must be freed to avoid memory leaks.
### SEE ALSO
[mm\_write\_banner()](#mm_write_banner),
 [mm\_read\_mtx\_crd\_size()](#mm_read_mtx_crd_size),
 [mm\_read\_mtx\_array\_size()](#mm_read_mtx_crd_size),
 [mm\_typecode\_to\_str()](#mm_typecode_to_str),
 [mm\_write\_mtx\_array\_size()](#mm_write_mtx_array_size),
---
The Matrix Market is a service of the
[Mathematical and Computational Sciences Division](http://math.nist.gov/mcsd/) /
[Information Technology Laboratory](http://www.itl.nist.gov/) /
[National Institute of Standards and Technology](http://www.nist.gov/)
[ [Home](http://math.nist.gov/MatrixMarket/index.html) ]
[ [Search](http://math.nist.gov/MatrixMarket/search.html) ]
[ [Browse](http://math.nist.gov/MatrixMarket/browse.html) ]
[ [Resources](http://math.nist.gov/MatrixMarket/resources.html) ]
Last change in this page : *May 2, 2000*.
[ <!--
enist('matrixmarket');
// -->
matrixmarket@nist.gov
].

