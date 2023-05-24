target_sources(
  micm_object
  PRIVATE
    interface.cpp
    interface.F90
)

target_include_directories(
  micm_object
  PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
)

target_include_directories(
  micm
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:.>
)

set_target_properties(micm_object PROPERTIES
  Fortran_MODULE_DIRECTORY ${MICM_MOD_DIR}
)