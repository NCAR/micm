################################################################################
# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)

# target_compile_features(micm INTERFACE cxx_std_20)

target_include_directories(micm
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(ENABLE_JSON)
  target_link_libraries(micm INTERFACE nlohmann_json::nlohmann_json)
  target_compile_definitions(micm INTERFACE USE_JSON)
endif()

if(ENABLE_OPENMP)
  target_link_libraries(micm INTERFACE OpenMP::OpenMP_CXX)
endif()

if(ENABLE_MPI)
  target_link_libraries(micm INTERFACE MPI::MPI_CXX)
endif()

if(ENABLE_CUDA OR ENABLE_OPENACC)
  add_library(micm_gpu)
  add_library(musica::micm_gpu ALIAS micm_gpu)

  target_link_libraries(micm_gpu 
    PRIVATE micm
  )

  if(ENABLE_CUDA)
    target_link_libraries(micm_gpu 
      PUBLIC cudart
    )
  endif()

  if(ENABLE_OPENACC)
    target_link_libraries(micm_gpu 
      PUBLIC OpenACC::OpenACC_CXX
    )
  endif()

  add_subdirectory(process)
endif()
