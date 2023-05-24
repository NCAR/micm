################################################################################
# object library
add_library(micm_object OBJECT)

target_include_directories(micm_object
  INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

# micm version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/micm/version.hpp @ONLY)

################################################################################
# micm solver
add_library(micm $<TARGET_OBJECTS:micm_object>)
add_library(musica::micm ALIAS micm)

if(ENABLE_JSON)
  target_link_libraries(micm INTERFACE nlohmann_json::nlohmann_json)
  target_compile_definitions(micm INTERFACE USE_JSON)
endif()

if(ENABLE_FORTRAN_LINK)
  set_target_properties(micm PROPERTIES
    Fortran_MODULE_DIRECTORY ${MICM_MOD_DIR}
  )

  target_include_directories(micm
    INTERFACE
      $<INSTALL_INTERFACE:include>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    PUBLIC  
      $<BUILD_INTERFACE:${MICM_MOD_DIR}>
      $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>
  )
  add_subdirectory(fortran_interface)
endif()

if(ENABLE_OPENMP)
  target_link_libraries(micm PRIVATE OpenMP::OpenMP_CXX)
endif()

if(ENABLE_MPI)
  target_link_libraries(micm PRIVATE MPI::MPI_CXX)
endif()
