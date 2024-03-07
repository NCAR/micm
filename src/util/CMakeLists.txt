if(MICM_ENABLE_CUDA)
  target_sources(micm_cuda
       PRIVATE
       cuda_vector_matrix.cu
  )
endif()
