include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
set(INSTALL_PREFIX "micm-${PROJECT_VERSION}" )

# micm_kokkos must ship in the same export set as micm: it is named in micm's
# interface by downstream projects, and an unexported target makes the installed
# package unusable for the Kokkos backend.
set(MICM_INSTALL_TARGETS micm)
if(MICM_ENABLE_KOKKOS)
  list(APPEND MICM_INSTALL_TARGETS micm_kokkos)
endif()

install(
  TARGETS
    ${MICM_INSTALL_TARGETS}
  EXPORT 
    micm_Exports
  LIBRARY DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  DIRECTORY
    ${PROJECT_SOURCE_DIR}/include/micm
  DESTINATION
    ${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
)

# Ship the license and third-party notices with the install tree.
install(
  FILES
    ${PROJECT_SOURCE_DIR}/LICENSE
    ${PROJECT_SOURCE_DIR}/NOTICE
  DESTINATION
    ${INSTALL_PREFIX}
)

# install the cmake config files
set(cmake_config_install_location ${INSTALL_PREFIX}/cmake)

install(
  EXPORT 
    micm_Exports
  DESTINATION 
    ${cmake_config_install_location}
  NAMESPACE musica::
)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/micmConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/micmConfig.cmake"
  INSTALL_DESTINATION
    ${cmake_config_install_location}
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/micmConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/micmConfig.cmake
    ${PROJECT_BINARY_DIR}/micmConfigVersion.cmake
  DESTINATION
    ${cmake_config_install_location}
)

######################################################################
# uninstall target

# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()