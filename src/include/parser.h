#ifndef PARSER
#define PARSER

#include "mmio.h" 
#include "sparseMatrix.h" 

typedef struct{
    uint row;
    uint col;
    double val;
} entry;     //MatrixMarket COO entry

/*
 * Parse MatrixMarket matrix stored in file at @matPath
 * Returns: allocated spmat sparse matrix with all field allocated
 * symetric matrix are expanded in a full matrix
 */
spmat* MMtoCSR(char* matPath);

/*
 * basic check for sparse matrix compliance to the app
 * posix bool return
 */
int MMCheck(MM_typecode typecode);
double* readVector(char* fpath);

#endif
