if(MICM_ENABLE_CUDA)
  target_sources(micm_cuda
     PRIVATE
     process_set.cu
  )
endif()
