#include <stdlib.h>

namespace micm
{

#ifdef __cplusplus
  extern "C"
  {
#endif

    typedef void (*FuncPtr)(double[], int64_t, int64_t);

    FuncPtr get_solver(char filepath[]);

#ifdef __cplusplus
  }  // extern "C"
#endif

}  // namespace micm