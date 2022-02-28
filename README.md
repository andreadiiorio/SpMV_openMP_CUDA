SpMV - Sparse parallel Matrix Vector Multiplication
========================================================
Andrea Di Iorio


Several implementations of Sparse parallel MatrixVector Multiplication in openMP and CUDA
#Implementations brief description
following an incremental numbering scheme (also used in .tex report and presentation)
the implementations are:
##openMP:
  CSR format implementations
    -sgemvSerial			serial implementations
    -spmvRowsBasicCSR		1 row per thread
    -spmvRowsBlocksCSR		1 block of row per thread
    -spmvTilesCSR			1 2D sub-block of the matrix per thread, inplace partitionated
    -spmvTilesAllocdCSR		1 2D sub-block of the matrix per thread, separate CSR per col partition
  ELL format implementations
    -spmvRowsBasicELL		1 row per thread
    -spmvRowsBlocksELL		1 row block per thread
    -spmvTilesELL			1 2D sub-block per thread
##CUDA
  CSR format implementations
	-cudaSpMVRowsCSR		1 thread per row
	-cudaSpMVWarpPerRowCSR  1 warp per row
  ELL format implementations
	-cudaSpMVRowsELL		1 thread per row, transposing the matrix and pitching for coalescing
	-cudaSpMVWarpsPerRowELLNTrasposed	1 warp per row (without trasposing since not necessary)

#Configurations evaluated
Beside the different implementations/partitioning schemes for the matrix
I've evaluated these additional configuarion, both 
##launch time 
-partitioning grid size (`gridRows`x`gridRows)
-scheduling configuration via omp env var:
	static
	dynamic (with chunk size adapatiation of static's chunk size / FAIR_CHUNKS_FOLDING(4)
##compile time
auxiliary rows lens vector (useful for premature ending of (parzial) point product accumulation with ELL's padding)
SIMD reduction (not usefull because of sparsity)
#Compilation
single main files builded separatelly with gcc for openMP only implementations
and with nvcc (unfortunatelly c++ frontend) for both openMP and CUDA implementations 
using __CUDACC__ nvcc exported macro to include in nvcc compilation, C extended statements for CUDA implementations
#Testing
threshold based confront of numerical result with a serial implementation, 
validating it with confront with CBLAS numerical result (previous a dense transformation of the problem)

See test/
#Perfomarce
TODO include pdfs in doc/ or tables in test/
