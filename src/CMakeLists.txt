######################################################################
# MICM library

set(MICM_SRC ODE_solver.F90
             ODE_solver_factory.F90
             ODE_solvers/mozart.F90
             ODE_solvers/rosenbrock.F90
             environment.F90
             kinetics.F90
             rate_constant_functions.F90
             preprocessor_output/factor_solve_utilities.F90
             preprocessor_output/kinetics_utilities.F90
             preprocessor_output/rate_constants_utility.F90
             core.F90)

add_library(micmlib ${MICM_SRC})

target_link_libraries(micmlib musicalib)

set_target_properties(micmlib PROPERTIES OUTPUT_NAME micm)
