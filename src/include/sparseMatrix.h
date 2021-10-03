//sparse matrix def & aux
//TODO adapt to work on both CUDA->ELL and std CSR
#ifndef SPARSEMATRIX
#define SPARSEMATRIX 
typedef struct{
    uint NZ,M,N;
    uint *JA;
    uint *IRP;
    //TODO MACRO EXPANSION TO FILTER CUDA STUFF 
    //uint MAXNZ;

    double *AS; 
} spmat; //describe a sparse matrix

#endif
