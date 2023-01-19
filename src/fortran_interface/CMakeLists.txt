target_sources(
  micm
  PRIVATE
    interface.cpp
)

target_include_directories(
  micm
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:.>
)