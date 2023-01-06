#include <micm/solver/chapman_ode_solver.hpp>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void solve(double* concentrations, uint64_t size);

#ifdef __cplusplus
}
#endif