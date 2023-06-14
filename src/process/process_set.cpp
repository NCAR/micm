#include <iostream>

namespace micm
{
  namespace openacc
  {

    void deriv()
    {
      const int N = 100000000;
      int data[N];

      // Initialize data
      for (int i = 0; i < N; ++i)
      {
        data[i] = i;
      }

// Parallelize the loop using OpenACC
#pragma acc parallel loop
      for (int i = 0; i < N; ++i)
      {
        data[i] *= 2;
      }
    }

  }  // namespace openacc
} // namespace micm