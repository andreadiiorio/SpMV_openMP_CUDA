#ifndef PARSER
#define PARSER

#include "mmio.h" 
#include "sparseMatrix.h" 

typedef struct{
    ulong row;
    ulong col;
    double val;
} entry;     //MatrixMarket COO entry


/* 
 * parse MatrixMarket matrix entries in @fp, of type @mcode and @NZ entries
 * into COOrdinate list of entries
 *  -> expand simmetric matrixes into a normal matrix with both parts
 *      so NZ will be inplace doubled
 * return allocated and filled COO entries with the NNZ number into
 */
entry* MMtoCOO(ulong* NZ, FILE *fp, MM_typecode mcode);
/*
 * write COO entries in @entries inside sparse matrix @mat
 * expected CSR arrays allocated
 * [simmetrical parts explicitly rappresented --> not important here]
 */
int COOtoCSR(entry* entries, spmat* mat);
/*
 * Parse MatrixMarket matrix stored in file at @matPath
 * IMPLEMENTED WRAPPING: MMtoCOO -> COOtoCSR
 * Returns: allocated spmat sparse matrix with all field allocated
 * symetric matrix are expanded in a full matrix
 */
spmat* MMtoCSR(char* matPath);

/*
 * basic check for sparse matrix compliance to the app
 * posix bool return
 */
int MMCheck(MM_typecode typecode);



#endif
