#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

char* SCHEDULES[]={"OMP_SCHED_STATIC","OMP_SCHED_DYNAMIC","OMP_SCHED_GUIDED","OMP_SCHED_AUTO"};
void ompGetRuntimeSchedule(int* kind_monotonic_chunk){
    /*
     * export OMP_SCHEDULE="[modifier:]kind[, chunk]"
     * typedef enum omp_sched_t {
     * //schedule kinds
     * omp_sched_static     = 0x1,
     * omp_sched_dynamic    = 0x2,
     * omp_sched_guided     = 0x3,
     * omp_sched_auto       = 0x4,
     * //schedule modifier
     * omp_sched_monotonic  = 0x80000000u
     * } omp_sched_t;
     */
    omp_sched_t k,kind; int chunk_size,monotonic;
    omp_get_schedule(&kind,&chunk_size);
    monotonic = omp_sched_monotonic & kind;
    k=kind;
    if(monotonic)   k = kind - omp_sched_monotonic;
    printf("omp sched gather:\tkind: %s - monotonic: %s\tomp chunkSize: %d\n",
      SCHEDULES[k-1],monotonic?"Y":"N",chunk_size);
    if (kind_monotonic_chunk){
        kind_monotonic_chunk[0] = k;
        kind_monotonic_chunk[1] = monotonic;
        kind_monotonic_chunk[2] = chunk_size;
    }
}

//WRAPPER to print all ICV vars
#ifdef OMP_GET_ICV_MAIN
int main(){
#else
void ompGetAllICV(){
#endif
    ompGetRuntimeSchedule(NULL);
}
