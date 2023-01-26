######################################################################
# MICM library

set(MICM_SRC constants.F90
  #             core.F90
             ODE_solver.F90
             ODE_solver_factory.F90
             ODE_solvers/mozart.F90
             ODE_solvers/rosenbrock.F90
             environment.F90
             kinetics.F90
             preprocessor_output/factor_solve_utilities.F90
             preprocessor_output/kinetics_utilities.F90
             preprocessor_output/rate_constants_utility.F90
             rate_constant.F90
             rate_constants/arrhenius.F90
             rate_constants/photolysis.F90
             rate_constants/ternary_chemical_activation.F90
             rate_constants/troe.F90
             rate_constants/wennberg_alkoxy.F90
             rate_constants/wennberg_nitrate.F90
             rate_constants/wennberg_tunneling.F90
             )

add_library(micmlib ${MICM_SRC})

target_link_libraries(micmlib musicalib)

set_target_properties(micmlib PROPERTIES OUTPUT_NAME micm)

######################################################################
# Set compiler options

if (ENABLE_OPENACC)

set(FLAGGED_FILES preprocessor_output/kinetics_utilities.F90
                  preprocessor_output/factor_solve_utilities.F90
                  preprocessor_output/rate_constants_utility.F90
                  rate_constants/arrhenius.F90
                  rate_constants/photolysis.F90
                  kinetics.F90
                  constants.F90   
                  ODE_solvers/rosenbrock.F90)

set_source_files_properties(${FLAGGED_FILES} PROPERTIES COMPILE_FLAGS "-acc -gpu=cc70,lineinfo,nofma,pinned -Minfo=all")

endif()
