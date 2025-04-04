if(MICM_ENABLE_CUDA)
  target_sources(micm_cuda
     PRIVATE
     lu_decomposition_mozart_in_place.cu
     linear_solver.cu
     rosenbrock.cu
  )
endif()