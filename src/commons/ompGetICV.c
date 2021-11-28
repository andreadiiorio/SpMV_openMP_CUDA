#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

char* SCHEDULES[]={"OMP_SCHED_STATIC","OMP_SCHED_DYNAMIC","OMP_SCHED_GUIDED","OMP_SCHED_AUTO"};
void ompGetRuntimeSchedule(int* kind_chunk_monotonic){
    /*
     * export OMP_SCHEDULE="[modifier:]kind[, chunk]"
     * typedef enum omp_sched_t {
     * //schedule kinds
     * omp_sched_static     = 0x1,
     * omp_sched_dynamic    = 0x2,
     * omp_sched_guided     = 0x3,
     * omp_sched_auto       = 0x4,
     * //schedule modifier  //TODO API>=5.0
     * omp_sched_monotonic  = 0x80000000u   //TODO OMP >=5
     * } omp_sched_t;
     */
    omp_sched_t k,kind; int chunk_size,monotonic=0;
    omp_get_schedule(&kind,&chunk_size);
    k=kind; //[monotonic OFF]
    #if _OPENMP >= 201811   //OMP_SCHED_SCHEDULE modifier from 5.0
    monotonic = omp_sched_monotonic & kind;
    if(monotonic)   k = kind - omp_sched_monotonic;
    #endif
    printf("omp sched gather:\tkind: %s\tomp chunkSize: %d\tmonotonic: %s\n",
      SCHEDULES[k-1],chunk_size,monotonic?"Y":"N");
    if (kind_chunk_monotonic){
        kind_chunk_monotonic[0] = k;
        kind_chunk_monotonic[1] = chunk_size;
        kind_chunk_monotonic[2] = monotonic;
    }
}

float ompVersionMacroMap(){
    switch ( _OPENMP ){
        case 200505:    return 2.5; 
        case 200805:    return 3.0;
        case 201107:    return 3.1;
        case 201307:    return 4.0;
        case 201511:    return 4.5;
        case 201811:    return 5.0;
        case 202011:    return 5.1;
    }
}
//WRAPPER to print all ICV vars
#ifdef OMP_GET_ICV_MAIN
int main(){
#else
void ompGetAllICV(){
#endif
    ompGetRuntimeSchedule(NULL);
    printf("omp API version:\t %1.1f\n",ompVersionMacroMap()); 
}
