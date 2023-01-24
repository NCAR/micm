target_sources(
  micm
  PRIVATE
    interface.cpp
    interface.F90
)

target_include_directories(
  micm
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:.>
)

set_target_properties(micm PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
)