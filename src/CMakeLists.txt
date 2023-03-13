################################################################################
# header only library
add_library(micm_header_only INTERFACE)
add_library(musica::micm_header_only ALIAS micm_header_only)

# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

target_include_directories(
  micm_header_only
  INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

################################################################################
# micm solver
add_library(micm)
add_library(musica::micm ALIAS micm)

target_link_libraries(
  micm
  PUBLIC
    micm_header_only
)

if(ENABLE_FORTRAN_LINK)
  add_subdirectory(fortran_interface)
endif()

if(ENABLE_OPENMP)
  target_link_libraries(micm PRIVATE OpenMP::OpenMP_CXX)
endif()

if(ENABLE_MPI)
  target_link_libraries(micm PRIVATE MPI::MPI_CXX)
endif()

if(ENABLE_JSON)
  target_link_libraries(micm PRIVATE nlohmann_json::nlohmann_json)
  target_compile_definitions(micm PUBLIC USE_JSON)
endif()