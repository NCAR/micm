################################################################################
# micm
add_library(micm)
add_library(musica::micm ALIAS micm)

# micm version
configure_file(version.cpp.in version.cpp @ONLY)

set_target_properties(
  micm 
  PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${MICM_LIB_DIR}
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

if(ENABLE_OPENMP)
  target_link_libraries(micm PUBLIC OpenMP::OpenMP_CXX)
endif()

if(ENABLE_MPI)
  target_link_libraries(micm PUBLIC MPI::MPI_CXX)
endif()

target_include_directories(
  micm
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

target_link_libraries(
  micm
  PRIVATE
)

target_sources(
  micm
  PRIVATE 
    ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)

add_subdirectory(system)