################################################################################
# micm
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)

# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

# if(ENABLE_OPENMP)
#   target_link_libraries(micm PUBLIC OpenMP::OpenMP_CXX)
# endif()

# if(ENABLE_MPI)
#   target_link_libraries(micm PUBLIC MPI::MPI_CXX)
# endif()

target_include_directories(
  micm
  INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)
