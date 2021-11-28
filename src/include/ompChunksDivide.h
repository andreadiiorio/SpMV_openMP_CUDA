#ifndef OMPCHUNKSDIVIDE
#define OMPCHUNKSDIVIDE

/*
 * chunks distribution for runtime schedule setted to dynamic, TODO guided
 * generic interface with elements, srcMatrix, configuration
 * TODO the srcMatrix is leaved for advanced chunk distribution...
 * configuration is expected to have a valid number of threadNum setted
 */
#include "config.h"
//distribution of @rows|blocks of @matrix, exploiting @config
typedef void (CHUNKS_DISTR )           (ulong,spmat*,CONFIG*);
typedef void (*CHUNKS_DISTR_INTERF )   (ulong,spmat*,CONFIG*);

//NOOP chunks division for manual configuration via export OMP_SCHEDULE  
inline void chunksNOOP(ulong r,spmat* mat,CONFIG* cfg){ return; }
#include "ompGetICV.h"
//fair division of @r elements from matrix @mat with threads num in @cfg
inline void chunksFair(ulong r,spmat* mat,CONFIG* cfg){
    assert(cfg->threadNum > 0); //configured target thread num
    omp_sched_t k,kind; int chunk_size,monotonic;
    omp_get_schedule(&kind,&chunk_size);
    monotonic = omp_sched_monotonic & kind;
    k = kind;
    if(monotonic)   k = kind-omp_sched_monotonic;
    switch(k){
        case omp_sched_static :
            DEBUG   printf("static it's already fair\n");
            return;
        case omp_sched_dynamic: 
            chunk_size = MAX(1,r/cfg->threadNum);
            break; 
        //case omp_sched_guided :
        //case omp_sched_auto   :
        //default:
    }
    omp_set_schedule(kind,chunk_size);
    DEBUG{  //check with ICV get
        printf("spmvRowsChunkFair:\tchunk adapted to %d",chunk_size);
        ompGetRuntimeSchedule(NULL);
    }
}
//fair division of @r elements from matrix @mat of threads in @cfg
//subdividing the fair share with a factor of @FAIR_CHUNKS_FOLDING
inline void chunksFairFolded(ulong r,spmat* mat,CONFIG* cfg){
    assert(cfg->threadNum > 0); //configured target thread num
    omp_sched_t k,kind; int chunk_size,monotonic;
    omp_get_schedule(&kind,&chunk_size);
    monotonic = omp_sched_monotonic & kind;
    k = kind;
    if(monotonic)   k = kind-omp_sched_monotonic;
    switch(k){
        case omp_sched_static :
            DEBUG   printf("static it's already fair\n");
            return;
        case omp_sched_dynamic: 
            chunk_size = MAX(1,r/(cfg->threadNum*FAIR_CHUNKS_FOLDING));
            break; 
        //case omp_sched_guided :
        //case omp_sched_auto   :
        //default:
    }
    omp_set_schedule(kind,chunk_size);
    DEBUG{  //check with ICV get
        printf("spmvRowsChunkFair:\tchunk adapted to %d",chunk_size);
        ompGetRuntimeSchedule(NULL);
    }
}
#endif
