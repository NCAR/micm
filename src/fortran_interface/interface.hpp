#include <stdlib.h>

namespace micm
{

#ifdef __cplusplus
  extern "C"
  {
#endif

    typedef void (*FuncPtr)(double[], uint64_t, uint64_t);

    FuncPtr get_solver(char filepath[]);

#ifdef __cplusplus
  }  // extern "C"
#endif

}  // namespace micm