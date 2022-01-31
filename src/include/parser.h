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

#ifndef PARSER
#define PARSER

#include "mmio.h" 
#include "sparseMatrix.h" 

typedef struct{
    ulong row;
    ulong col;
    double val;
} entry;     //MatrixMarket COO entry

typedef struct{
    MM_typecode mcode;
    entry* entries;
    ulong* rowLens;
    ulong M,N,NZ;   //spmat sizes
} MatrixMarket;

////COO PARSE
//parse and check MatrixMarket mat in @matPath file
MatrixMarket* MMRead(char* matPath);
void freeMatrixMarket(MatrixMarket* mm);
//basic check for sparse matrix compliance to the app, return posix bool
int MMCheck(MM_typecode typecode);
/* 
 * parse MatrixMarket matrix entries in @fp, of type @mcode
 * into COOrdinate list of entries
 *  -> expand simmetric matrixes into a normal matrix with both parts
 *      so NZ will be inplace doubled
 * return allocated and filled COO entries with the NNZ number into
 * NO SORT CHECKING HERE
 */
entry* MMtoCOO(ulong* NZ, FILE *fp, MM_typecode mcode,ulong* rowLens);

////COO -> ANYTHING ELSE CONVERSION
/*
 * write COO entries in @entries inside sparse matrix @mat in ELL format
 * EXPECTED: CSR arrays allocated, @entries col sorted in (not madatory consecut) rows
 * [simmetrical parts explicitly rappresented --> not important here]
 */
int COOtoCSR(entry* entries, spmat* mat,ulong* rowLens);
/*
 * write COO entries in @entries inside sparse matrix @mat in ELL format
 * EXPECTED: @entries col sorted in (not madatory consecut) rows
 * ELL internal array allocated in this function, not freed in case of error
 */
int COOtoELL(entry* entries, spmat* mat, ulong* rowLens);
////wrapper MM -> specialized target
/*
 * Parse MatrixMarket matrix stored in file at @matPath
 * IMPLEMENTED WRAPPING: MMtoCOO -> COOtoCSR
 * Returns: allocated spmat sparse matrix with all field allocated
 * symetric matrix are expanded in a full matrix
 */
spmat* MMtoCSR(char* matPath);
spmat* MMtoELL(char* matPath);


#endif
