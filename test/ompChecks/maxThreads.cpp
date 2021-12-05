#include <iostream>
#include <string>
#include <vector>
#include "omp.h"
int main() {

omp_set_num_threads(1);

std::cout << "before parallel section: " << std::endl;
std::cout << "Num threads = " << omp_get_num_threads() << std::endl;
std::cout << "Max threads = " << omp_get_max_threads() << std::endl;

//without num_threads(5), only 1 thread is created
#pragma omp parallel num_threads(5) 
  {          
      #pragma omp single
      {
          std::cout << "inside parallel section: " << std::endl;
          std::cout << "Num threads = " << omp_get_num_threads() << std::endl;
          std::cout << "Max threads = " << omp_get_max_threads() << std::endl;
      }
  }

  return 0;
}
