#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(){
    omp_sched_t kind; int chunk_size;
    omp_get_schedule(&kind,&chunk_size);
    printf("omp sched kind:\t%d omp chunkSize:\t%d\n",kind,chunk_size);
    if (kind == omp_sched_static)   printf("STATIC");
}
