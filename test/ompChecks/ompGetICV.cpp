#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <iostream>
    
int main(){
    printf("export OMP_DISPLAY_ENV=VERBOSE for full ICV details\n"); 
    printf("omp MAX THREADS USABLE\t%d\n",omp_get_max_threads());
    std::cout << "In parallel? " << omp_in_parallel() << std::endl;
    std::cout << "Num threads = " << omp_get_num_threads() << std::endl;
    std::cout << "Max threads = " << omp_get_max_threads() << std::endl;
    
    std::cout << "Entering my parallel region: " << std::endl;
    
    //without num_threads(5), only 1 thread is created
    #pragma omp parallel num_threads(5) 
      {          
          #pragma omp single nowait
          {
              std::cout << "In parallel? " << omp_in_parallel() << std::endl;
              std::cout << "Num threads = " << omp_get_num_threads() << std::endl;
              std::cout << "Max threads = " << omp_get_max_threads() << std::endl;
          }
      }
}
