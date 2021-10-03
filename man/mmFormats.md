Text File Formats
=================
We briefly describe the ASCII file formats for matrices redistributed by the Matrix Market :
- [Matrix Market Exchange Formats](#MMformat)
- [Harwell-Boeing Exchange Format](#hb)
- [Coordinate Text File Format](#coord) (**to be phased out**)

#Matrix Market Exchange Formats
A complete description is provided in the paper  *The Matrix Market Formats: Initial Design* reports/MMformat.ps.gz

MM exchange format -> Collection of affiliated formats which share design elements. Initial specification of 2 formats.
##Coordinate Format
For representing general sparse matrices.
Only nonzero entries are provided, and the coordinates of each nonzero entry is given explicitly. 

example of a real 5x5 general sparse matrix.
```

             1    0      0       6      0     
             0   10.5    0       0      0     
             0    0    .015      0      0     
             0  250.5    0     -280    33.32  
             0    0      0       0     12     
```

In MM **Coordinate Format** this could be represented as follows.
```

%%MatrixMarket matrix coordinate real general
%=================================================================================
%
% This ASCII file represents a sparse MxN matrix with L 
% nonzeros in the following Matrix Market format:
%
% +----------------------------------------------+
% |%%MatrixMarket matrix coordinate real general | <--- header line
% |%                                             | <--+
% |% comments                                    |    |-- 0 or more comment lines
% |%                                             | <--+         
% |    M  N  L                                   | <--- rows, columns, entries
% |    I1  J1  A(I1, J1)                         | <--+
% |    I2  J2  A(I2, J2)                         |    |
% |    I3  J3  A(I3, J3)                         |    |-- L lines
% |        . . .                                 |    |
% |    IL JL  A(IL, JL)                          | <--+
% +----------------------------------------------+   
%
% Indices are 1-based, i.e. A(1,1) is the first element.
%
%=================================================================================
  5  5  8
    1     1   1.000e+00
    2     2   1.050e+01
    3     3   1.500e-02
    1     4   6.000e+00
    4     2   2.505e+02
    4     4  -2.800e+02
    4     5   3.332e+01
    5     5   1.200e+01
```

The first line contains the type code. In this example, it indicates that the
object being represented is a matrix in coordinate format and that the numeric
data following is real and represented in general form. (By general we mean
that the matrix format is not taking advantage of any symmetry properties.)

Variants of the coordinate format are defined for matrices with complex and
integer entries, as well as for those in which only the position of the nonzero
entries is prescribed (pattern matrices). (These would be indicated by
changing real to complex, integer, or
pattern, respectively, on the header line). 
###pattern smart variants are
defined for cases in which **symmetries** can be used to significantly reduce the
size of the data: *symmetric, skew-symmetric* and *Hermitian*. In these cases, **only
entries in the lower triangular portion need be supplied**. In the skew-symmetric
case the diagonal entries are zero, and hence they too are omitted. (These
would be indicated by changing general to symmetric,
skew-symmetric, or hermitian, respectively, on the header
line).
    
###Array Format
    For representing general dense matrices.
    All entries are provided in a pre-defined (column-oriented) order.

#Harwell-Boeing Exchange Format
------------------------------

Most popular mechanism for text-file exchange of sparse matrix data.
fullSpec:   [User's Guide for the Harwell-Boeing Sparse Matrix Collection](ftp://ftp.cerfacs.fr/pub/algo/matrices/harwell_boeing/userguide.ps.Z)


- Matrix data is held in an 80-column, fixed-length format for portability. 
- Multiple line header block, which is followed by two, three, or four data blocks. 
    The header block contains summary information on the storage formats and space requirements. 

If there are no right-hand-side vectors, the matrix has
a four-line header block followed by two or three data blocks
containing, in order, the column (or element) start
pointers, the row (or variable) indices, and the numerical
values. If right-hand sides are present, there is a fifth
line in the header block and a fourth data block
containing the right-hand side(s). The blocks containing
the numerical values and right-hand side(s) are optional.
The right-hand side(s) can be present only when the
numerical values are present.
If right-hand sides are present, then vectors for starting guesses
and the solution can also be present; if so, they appear as separate
full arrays in the right-hand side block following the right-hand
side vector(s).

The first line contains the 72-character title and the
8-character identifier by which the matrix is referenced
in our documentation.
The second line contains the number of lines for each of the
following data blocks as well as the total number of lines,
excluding the header block. The third line
contains a three character string denoting the matrix type
as well as the number of rows, columns (or elements),
entries, and, in the case of unassembled matrices, the total
number of entries in elemental matrices. The
fourth line contains the variable Fortran formats for the
following data blocks. The fifth line is
present only if there are right-hand sides. It contains a one
character string denoting the storage format for the
right-hand sides as well as the number of right-hand sides,
and the number of row index entries (for the assembled case).
The exact format is given by the following, where the names of the
Fortran variables in the subsequent programs are given in parenthesis:
  
  
**Line 1** `(A72,A8)`
|  |  |
| --- | --- |
| 
Col. 1 - 72 
 | 
Title (`TITLE`)
 |
| 
Col. 73 - 80 
 
Key (`KEY`)
 | |
  
**Line 2** `(5I14)`
|  |  |
| --- | --- |
| 
Col. 1 - 14 
 | 
Total number of lines excluding header (`TOTCRD`)
 |
| 
Col. 15 - 28 
 | 
Number of lines for pointers (`PTRCRD`)
 |
| 
Col. 29 - 42 
 | 
Number of lines for row (or variable) indices (`INDCRD`)
 |
| 
Col. 43 - 56 
 | 
Number of lines for numerical values (`VALCRD`)
 |
| 
Col. 57 - 70 
 | 
Number of lines for right-hand sides (`RHSCRD`)
 |
|  | 
(including starting guesses and solution vectors if present)
 |
|  | 
(zero indicates no right-hand side data is present)
 |
  
 **Line 3** `(A3, 11X, 4I14)`
|  |  |
| --- | --- |
| 
Col. 1 - 3 
 | 
Matrix type (see below) (`MXTYPE`)
 |
| 
Col. 15 - 28 
 | 
Number of rows (or variables) (`NROW`)
 |
| 
Col. 29 - 42 
 | 
Number of columns (or elements) (`NCOL`)
 |
| 
Col. 43 - 56 
 | 
Number of row (or variable) indices (`NNZERO`)
 |
|  | 
(equal to number of entries for assembled matrices)
 |
|  |  |
| 
Col. 57 - 70 
 | 
Number of elemental matrix entries (`NELTVL`)
 |
|  | 
(zero in the case of assembled matrices)
 |
  
**Line 4** `(2A16, 2A20)`
|  |  |
| --- | --- |
| 
Col. 1 - 16 
 | 
Format for pointers (`PTRFMT`)
 |
| 
Col. 17 - 32 
 | 
Format for row (or variable) indices (`INDFMT`)
 |
| 
Col. 33 - 52 
 
Format for numerical values of coefficient matrix (`VALFMT`)
 | |
| 
Col. 53 - 72 
 | 
Format for numerical values of right-hand sides (`RHSFMT`)
 |
  
**Line 5** `(A3, 11X, 2I14)` *Only present if there are right-hand sides present*
|  |  |
| --- | --- |
| 
Col. 1 
 | 
Right-hand side type:
 |
|  | `F` for full storage or
 |
|  | `M` for same format as matrix
 |
| 
Col. 2 
 | `G` if a starting vector(s) (Guess) is supplied. (`RHSTYP`)
 |
| 
Col. 3 
 | `X` if an exact solution vector(s) is supplied.
 |
| 
Col. 15 - 28 
 | 
Number of right-hand sides (`NRHS`)
 |
| 
Col. 29 - 42 
 | 
Number of row indices (`NRHSIX`)
 |
|  | 
(ignored in case of unassembled matrices)
 |
**Note:**
*For matrices in elemental form, the leading two dimensions in
the header give the number of variables in the finite element
application and the number of elements. It is common that not all of
the variables in the application appear in the linear algebra
subproblem; hence the matrix represented can be of lower order than
the first parameter, described as the "number of variables*
(`NROW`)".
*The finite element variables are numbered from 1 to*
`NROW`, 
*but only the subset of variables that actually appear in the list of
variables for the elements define the rows and columns of the
matrix. The actual order of the square matrix cannot be determined
until all of the indices are read.*

The three character type field on line 3 describes the matrix type.
The following table lists the permitted values for each of the three
characters. As an example of the type field, RSA denotes that
the matrix is real, symmetric, and assembled.

  
  
**First Character:**
|  |
| --- |
| `R` Real matrix
 |
| `C` Complex matrix
 |
| `P` Pattern only (no numerical values supplied)
 |
  
**Second Character:**
|  |
| --- |
| `S` Symmetric
 |
| `U` Unsymmetric
 |
| `H` Hermitian
 |
| `Z` Skew symmetric
 |
| `R` Rectangular
 |
  
**Third Character:**
|  |
| --- |
| `A` Assembled
 |
| `E` Elemental matrices (unassembled)
 |
### Example Fortran Code for Reading Harwell-Boeing Files

To formalize the logical block structure of the data, we
have included two pieces of sample FORTRAN code for reading
a matrix in the format of the sparse matrix test collection.
Both codes
assume the data comes from input unit LUNIT.
Neither is a complete code. Real code should include error
checking to ensure that the target arrays into which the
data are read are large enough. The design allows the
arrays to be read by a separate subroutine that can avoid
the use of possibly inefficient implicit DO-loops.

* [First sample fortran code](src/hbcode1.f) :
the standard case, a sparse matrix in standard format with no right-hand sides.
* [Second sample fortran code](src/hbcode2.f) :
illustrates the full generality of the representation.

 The code above outlines the structure of the data. The
interpretation of the row (or variable) index arrays
will require knowledge of the matrix and right-hand side
types, as read in this code.

### Matlab Procedures for Reading/Writing Harwell-Boeing Files

The developers of the NEP matrix collection have provided a Matlab m-file to
write a [Matlab sparse matrix in Harwell-Boeing format](src/dm2hb.m).
A [version for complex matrices](src/zm2hb.m) is also available.

The Berkeley Benchmarking and Optimization (BeBOP) Group has developed a 
[library and standalone utility](http://www.nist.gov/cgi-bin/exit_nist.cgi?timeout=5&url=http://www.cs.berkeley.edu/~mhoemmen/bebop/smc.html)
for converting between Harwell-Boeing, Matrix Market,
and MATLAB sparse matrix formats.


#Coordinate Text File
--------------------
**Note:** __This format is being phased out!!!!!!!!!!!!__

- 1° linea: the number of rows *m*, columns *n*, and nonzeros *nz* in the matrix. 
- The nonzero matrix elements are then listed, one per line, 
    specifying row index *i*, column index *j*, and the value *a(i,j)*, in that order.
For example,
```

      m       m       nz
     i1      j1      val1
     i2      j2      val2
     i3      j3      val3
     .       .        .
     .       .        .
     .       .        .
     inz     jnz     valnz
```

White space is not significant, (i.e. a fixed column is not used). The nonzero
values may be in either in fixed or floating point representation, to any
precision (although Fortran and C typically parse less than 20 significant digits).
For example, the following are each acceptable: 3, 3.141,
+3.1415626536E000, 3.1e0.

Experiments show that these coordinate files are approximately 30% larger than
corresponding Harwell-Boeing files. Versions compressed with Unix compress 
or gzip typically exhibits similar ratios. 

To represent only structure information of a sparse matrix, a single zero can
be placed in the *value* position, e.g.

```

      M       N       nz
     i1      j1       0
     i2      j2       0
     i3      j3       0
     .       .        .
     .       .        .
     .       .        .
     inz     jnz      0
```


Although more efficient schemes are available, this allows the same routine to
read both types of files. The addition of a single byte to each line of the
file is typically of little consequence.

Note that there is no implied order for the matrix elements. This allows one
to write simple print routines which traverse the sparse matrix in whatever
natural order given by the particular storage scheme.

Also note that no annotations are used for storing matrices with special
structure. (This keeps the parsing routines simple.) Symmetric matrices can be
represented by only their upper or lower triangular portions, but the file
format reveals just that --- the reading program sees only a triangular
matrix. (The application is responsible for reinterpreting this.)
A [MATLAB](http://www.mathworks.com/) function ([M-file](src/rdcoord.m)) is available which reads a matrix in coordinate text file format and creates a sparse matrix is available.
---


#SW
- [Matrix Market I/O in C](mmio-c.html)
- [Matrix Market I/O in Fortran](mmio/f/mmiof77.html)
- [Matrix Market I/O in Matlab®](mmio/matlab/mmiomatlab.html)
- [BeBOP Sparse Matrix Conversion Library](http://www.nist.gov/cgi-bin/exit_nist.cgi?timeout=5&url=http://www.cs.berkeley.edu/~mhoemmen/bebop/smc.html)
- [Matrix Market I/O in Python](http://www.nist.gov/cgi-bin/exit_nist.cgi?timeout=5&url=http://docs.scipy.org/doc/scipy/reference/tutorial/io.html#matrix-market-files)
- [Matrix Market I/O in Gensim Python framework](http://www.nist.gov/cgi-bin/exit_nist.cgi?timeout=5&url=http://nlp.fi.muni.cz/projekty/gensim/tut1.html) (enables memory-efficient conversion to several other sparse formats)

