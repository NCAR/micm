#include <stdlib.h>

// namespace micm {

#ifdef __cplusplus
extern "C" {
#endif

  typedef void (*FuncPtr)(double [], double [], double []);

  FuncPtr get_solver(char filepath[]);

#ifdef __cplusplus
} // extern "C"
#endif

// }