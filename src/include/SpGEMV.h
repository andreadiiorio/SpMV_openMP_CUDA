#ifndef _SPGEMV
#define _SPGEMV

#include "macros.h"
#include "sparseMatrix.h"

#define _ROWS            "ROWS"                                                                                                                                                                                                                                               
#define _SORTED_ROWS     "SORTED_ROWS"                                                                                                                                                                                                                                        
#define _TILES           "TILES"                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                              
typedef enum {                                                                                                                                                                                                                                                                
    ROWS,                                                                                                                                                                                                                                                                     
    SORTED_ROWS,                                                                                                                                                                                                                                                              
    TILES                                                                                                                                                                                                                                                                     
} COMPUTE_MODE;                                                                                                                                                                                                                                                               


/*
 * basic spgemv with row partitioning, 1 row per thread in consecutive order 
 * statically assigned to threads
 */                                                                                                                                                                                                                                                                           
int spgemvRowsBasic(spmat* mat, double* vect, CONFIG* cfg, double* outVect);

//serial implementaion of sparse GEMV
int sgemvSerial(spmat* mat,double* vect, CONFIG* cfg, double* outVect);

typedef int (SPGEMV)         (spmat*,double*,CONFIG*,double*);
typedef int (*SPGEMV_INTERF) (spmat*,double*,CONFIG*,double*);

static const SPGEMV_INTERF  SpgemvFuncs[] = {
    &spgemvRowsBasic,
    NULL,
};
#endif
