#include "cudaUtils.h"
//int cudaErr(cudaError_t e,const char* errStr);
int spMatCpyCSR(spmat* m,spmat* dst){
	char* errS = "spMatCpyCSR";
	//prepare the dst-copy struct
	spmat dstLocal;
	memset(&dstLocal,0,sizeof(dstLocal));
	dstLocal.M  = m->M;
	dstLocal.N  = m->N;
	dstLocal.NZ = m->NZ;

	ulong nz = m->NZ;
	//allocs and copies for the CSR fields inside @m
	if(cudaErr( cudaMalloc( &dstLocal.JA,sizeof(*(dstLocal.JA))*nz), errS))								goto err;
	if(cudaErr( cudaMemcpy( dstLocal.JA,m->JA,sizeof(*(dstLocal.JA))*nz,dirUp), errS))					goto err;

	if(cudaErr( cudaMalloc( &dstLocal.AS,sizeof(*(dstLocal.AS)) * nz ), errS))							goto err;
	if(cudaErr( cudaMemcpy( dstLocal.AS,m->AS,sizeof(*(dstLocal.AS))*nz,dirUp), errS))					goto err;
	
	if(cudaErr( cudaMalloc( &dstLocal.IRP,sizeof(*(dstLocal.IRP)) * (dstLocal.M+1) ), errS))		  	goto err;
	if(cudaErr( cudaMemcpy( dstLocal.IRP,m->IRP,sizeof(*(dstLocal.IRP)) * (dstLocal.M+1),dirUp), errS)) goto err;
	#ifdef ROWLENS
	if(cudaErr( cudaMalloc( &dstLocal.RL,sizeof(*(dstLocal.RL)) * dstLocal.M ), errS))					goto err;
	if(cudaErr( cudaMemcpy( dstLocal.RL,m->RL,sizeof(*(dstLocal.RL)) * dstLocal.M,dirUp), errS))		goto err;
	#endif
	//write dst-copy in destination CUDA mem struct
	if(cudaErr( cudaMemcpy( dst,&dstLocal,sizeof(*dst),dirUp), errS))	goto err;
	return EXIT_SUCCESS;

	err:
	cudaFree(dstLocal.JA);
	cudaFree(dstLocal.AS);
	cudaFree(dstLocal.IRP);
	#ifdef ROWLENS
	cudaFree(dstLocal.RL);
	#endif
	return EXIT_FAILURE;
}
int spMatCpyELL(spmat* m,spmat* dst){
	char* errS = "spMatCpyELL";
	//prepare the dst-copy struct
	spmat dstLocal;
	memset(&dstLocal,0,sizeof(dstLocal));
	dstLocal.M = m->M;
	dstLocal.N = m->N;
	dstLocal.NZ= m->NZ;
	dstLocal.MAX_ROW_NZ= m->MAX_ROW_NZ;
	size_t pitchJAB,pitchASB;	//2D allocation rows sizes in bytes
	ulong nz = m->NZ;
	ulong maxRow = m->MAX_ROW_NZ;
	//sizes in Bytes
	size_t jaB = sizeof(m->JA), asB = sizeof(m->AS);
	//allocs and copies for the ELL fields inside @m
	if(cudaErr( cudaMallocPitch(&dstLocal.JA,&pitchJAB,maxRow*jaB,m->M), errS))						goto err;
	if(cudaErr( cudaMemcpy2D(dstLocal.JA,pitchJAB,m->JA,maxRow*jaB,maxRow*jaB,m->M,dirUp),errS))	goto err;
	
  	if(cudaErr( cudaMallocPitch(&dstLocal.AS,&pitchASB,maxRow*asB,m->M), errS))						goto err;
	if(cudaErr( cudaMemcpy2D(dstLocal.AS,pitchASB,m->AS,maxRow*asB,maxRow*asB,m->M,dirUp), errS))	goto err;
	#ifdef ROWLENS
	if(cudaErr( cudaMalloc(&dstLocal.RL,sizeof(*(dstLocal.RL)*dstLocal.M)), errS))					goto err;
	if(cudaErr( cudaMemcpy(dstLocal.RL,m->RL,sizeof(*(dstLocal.RL)*dstLocal.M)), errS))				goto err;
	#endif
	//set pitches as num of elements
	dstLocal.pitchJA = pitchJAB / sizeof(*(dstLocal.JA));
	dstLocal.pitchAS = pitchASB / sizeof(*(dstLocal.AS));
	DEBUG  	
		printf("spMatCpyELL\tcudaMallocPitch:\t\tAS: pitched %u extra elements\tJA: pitched %u extra elements\n",
			dstLocal.pitchAS - dstLocal.MAX_ROW_NZ,dstLocal.pitchJA - dstLocal.MAX_ROW_NZ);
	//write dst-copy in destination CUDA mem struct
	if(cudaErr( cudaMemcpy( dst,&dstLocal,sizeof(dstLocal),dirUp), errS))								goto err;
	return EXIT_SUCCESS;

	err:
	cudaFree(dstLocal.JA);
	cudaFree(dstLocal.AS);
	#ifdef ROWLENS
	cudaFree(dstLocal.RL);
	#endif
	return EXIT_FAILURE;
}

int spMatCpyELLNNPitched(spmat* m,spmat* dst){
	char* errS = "spMatCpyELLNNPitched";
	//prepare the dst-copy struct
	spmat dstLocal;
	memset(&dstLocal,0,sizeof(dstLocal));
	dstLocal.M = m->M;
	dstLocal.N = m->N;
	dstLocal.NZ= m->NZ;
	dstLocal.MAX_ROW_NZ= m->MAX_ROW_NZ;
	ulong nz = m->NZ;
	ulong maxRow = m->MAX_ROW_NZ;
	//sizes in Bytes
	size_t jaB = sizeof(m->JA), asB = sizeof(m->AS);
	//allocs and copies for the ELL fields inside @m
	if(cudaErr( cudaMalloc(&dstLocal.JA,maxRow*jaB*m->M), errS))						goto err;
	if(cudaErr( cudaMemcpy(dstLocal.JA,m->JA,maxRow*jaB*m->M,dirUp),errS))	goto err;
	
  	if(cudaErr( cudaMalloc(&dstLocal.AS,maxRow*asB*m->M), errS))						goto err;
	if(cudaErr( cudaMemcpy(dstLocal.AS,m->AS,maxRow*asB*m->M,dirUp), errS))	goto err;
	#ifdef ROWLENS
	if(cudaErr( cudaMalloc(&dstLocal.RL,sizeof(*(dstLocal.RL)*dstLocal.M)), errS))					goto err;
	if(cudaErr( cudaMemcpy(dstLocal.RL,m->RL,sizeof(*(dstLocal.RL)*dstLocal.M)), errS))				goto err;
	#endif
	//set pitches as num of elements
	dstLocal.pitchJA = dstLocal.MAX_ROW_NZ;
	dstLocal.pitchAS = dstLocal.MAX_ROW_NZ;
	DEBUG  	
		printf("spMatCpyELLNNPitched\tcudaMallocPitch:\t\tAS: pitched %u extra elements\tJA: pitched %u extra elements\n",
			dstLocal.pitchAS - dstLocal.MAX_ROW_NZ,dstLocal.pitchJA - dstLocal.MAX_ROW_NZ);
	//write dst-copy in destination CUDA mem struct
	if(cudaErr( cudaMemcpy( dst,&dstLocal,sizeof(dstLocal),dirUp), errS))								goto err;
	return EXIT_SUCCESS;

	err:
	cudaFree(dstLocal.JA);
	cudaFree(dstLocal.AS);
	#ifdef ROWLENS
	cudaFree(dstLocal.RL);
	#endif
	return EXIT_FAILURE;
}
