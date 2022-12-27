# micm version
configure_file(version.cpp.in version.cpp @ONLY)
add_library(micm_version ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)

target_include_directories(
  micm_version
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

set_target_properties(
  micm_version 
  PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${MICM_LIB_DIR}
)

# micm
add_library(micm)
add_library(musica::micm ALIAS micm)

set_target_properties(
  micm 
  PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${MICM_LIB_DIR}
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(
  micm
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

target_link_libraries(
  micm
  PRIVATE
    micm_version
)

target_sources(
  micm
  PRIVATE 
    dummy.cpp
)
