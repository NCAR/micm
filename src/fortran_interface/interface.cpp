#ifdef __cplusplus
extern "C" {
#endif

/// @brief A function that allows Fortran to get access to any micm solver
/// @param filepath A configuration file that defines the chemical system being solved
extern void get_solver(char filepath[]);

#ifdef __cplusplus
}
#endif