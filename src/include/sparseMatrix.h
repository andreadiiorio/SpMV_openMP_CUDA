//sparse matrix def & aux
//TODO adapt to work on both CUDA->ELL and std CSR
#ifndef SPARSEMATRIX
#define SPARSEMATRIX 
typedef struct{
    uint NZ,M,N;
    uint* JA;
    //CSR SPECIFIC
#ifdef ROWLENS
    uint* RL;   //row lengths
#endif
    uint* IRP;
    //CUDA SPECIFIC
    //uint MAXNZ;

    double *AS; 
} spmat; //describe a sparse matrix

#endif
