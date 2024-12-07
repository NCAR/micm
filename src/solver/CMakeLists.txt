if(MICM_ENABLE_CUDA)
  target_sources(micm_cuda
     PRIVATE
     lu_decomposition_doolittle.cu
     linear_solver.cu
     rosenbrock.cu
  )
endif()