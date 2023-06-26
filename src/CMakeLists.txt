################################################################################
# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm INTERFACE)
add_library(musica::micm ALIAS micm)
target_compile_features(micm INTERFACE cxx_std_20)

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
  target_link_libraries(micm PRIVATE OpenMP::OpenMP_CXX)
endif()

if(ENABLE_MPI)
  target_link_libraries(micm PRIVATE MPI::MPI_CXX)
endif()
