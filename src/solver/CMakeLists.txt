target_sources(micm_cuda
    PRIVATE
    lu_decomposition_mozart_in_place.cu
    lu_decomposition_mozart_in_place_hardcoded.cu
    linear_solver_in_place.cu
    rosenbrock.cu
)