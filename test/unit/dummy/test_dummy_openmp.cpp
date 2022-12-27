#include <micm/dummy.hpp>
#include <assert.h>
#include <omp.h>
#include <iostream>

int main(){
  const int nthreads = 4;
  bool failed[nthreads] = {false};

  omp_set_num_threads(nthreads);
  #pragma omp parallel
  {
    try{
        dummy();
    }
    catch(...)
    {
      failed[omp_get_thread_num()] = true;
    }
  }

  #pragma omp barrier

  for(uint8_t i = 0; i < nthreads; ++i)
  {
    assert(!failed[i]);
  }
}