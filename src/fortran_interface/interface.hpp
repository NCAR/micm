#include <stdlib.h>

// namespace micm {

  typedef void (*FuncPtr)(double *, double *, double *);

  FuncPtr get_solver(char filepath[]);
// }