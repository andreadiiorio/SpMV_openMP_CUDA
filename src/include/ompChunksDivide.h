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
    omp_sched_t k,kind; int chunk_size,chunk_size_new=0,monotonic;
    omp_get_schedule(&kind,&chunk_size);
    k = kind;
    #if _OPENMP >= 201811   //OMP_SCHED_SCHEDULE modifier from 5.0
    monotonic = omp_sched_monotonic & kind;
    if(monotonic)   k = kind-omp_sched_monotonic;
    #endif
    (void) monotonic;   //else no unused warning
    switch(k){
        case omp_sched_static :
            DEBUG   printf("static it's already fair\n");
            return;
        case omp_sched_dynamic: 
            chunk_size_new = MAX(1,r/cfg->threadNum);
            break; 
        //case omp_sched_guided :
        //case omp_sched_auto   :
        //default:
    }
    if(chunk_size == chunk_size_new)    return;
    omp_set_schedule(kind,chunk_size_new);
    VERBOSE printf("chunksFair:\tchunk adapted to %d\n",chunk_size_new);
    DEBUG   ompGetRuntimeSchedule(NULL);
}
//fair division of @r elements from matrix @mat of threads in @cfg
//subdividing the fair share with a factor of @FAIR_CHUNKS_FOLDING
inline void chunksFairFolded(ulong r,spmat* mat,CONFIG* cfg){
    assert(cfg->threadNum > 0); //configured target thread num
    omp_sched_t k,kind; int chunk_size,chunk_size_new=0,monotonic;
    omp_get_schedule(&kind,&chunk_size);
    k = kind;
    #if _OPENMP >= 201811   //OMP_SCHED_SCHEDULE modifier from 5.0
    monotonic = omp_sched_monotonic & kind;
    if(monotonic)   k = kind-omp_sched_monotonic;
    #endif
    (void) monotonic;   //else no unused warning
    switch(k){
        case omp_sched_static :
            DEBUG   printf("static it's already fair\n");
            return;
        case omp_sched_dynamic: 
            chunk_size_new = MAX(1,r/(cfg->threadNum*FAIR_CHUNKS_FOLDING));
            break; 
        //case omp_sched_guided :
        //case omp_sched_auto   :
        //default:
    }
    if(chunk_size == chunk_size_new)    return;
    omp_set_schedule(kind,chunk_size_new);
    DEBUG{  //check with ICV get
        printf("chunksFairFolded:\tchunk adapted to %d\n",chunk_size_new);
        ompGetRuntimeSchedule(NULL);
    }
}
#endif
