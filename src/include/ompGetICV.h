#ifndef OMPGETICV_H
#define OMPGETICV_H
//only header definitions
/*
 * log sched configuration on stdout
 * return kind,monotonic,chunkSize if arg0 not NULL
 */
void ompGetRuntimeSchedule(int* );
void ompGetAllICV();    //only if not OMP_GET_ICV_MAIN 
#endif
