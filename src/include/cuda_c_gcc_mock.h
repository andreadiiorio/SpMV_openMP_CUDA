/*
 * mock NVCC.C++ FRONTEND NEEDED STUFF FOR A GCC COMPILATION
 */
#ifndef CUDA_C_GCC_MOCK_H
#define CUDA_C_GCC_MOCK_H

#pragma message("compiling forwarded to HOST COMPILER")


////TODO MOCK NVCC C++ FRONTEND
typedef struct{
	double x;
	double y;
	double z;
} dim3;	///TODO USELESS??
//#define decltype(type)	void*	//TODO REMOVE? MOCK C++ DYN CAST 
#define decltype(type)	typeof(type)
#define cudaFree(x)		(void) 1

#endif
