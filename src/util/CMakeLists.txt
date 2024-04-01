if(MICM_ENABLE_CUDA)
  target_sources(micm_cuda
       PRIVATE
       cuda_matrix.cu
  )
endif()
