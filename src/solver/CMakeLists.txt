if(MICM_ENABLE_CUDA)
  target_sources(micm_cuda
     PRIVATE
     lu_decomposition.cu
     linear_solver.cu
     rosenbrock.cu
  )
endif()