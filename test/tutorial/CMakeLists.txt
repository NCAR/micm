################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME README_example SOURCES test_README_example.cpp)
create_standard_test(NAME rate_constants_no_user_defined_example_by_hand SOURCES test_rate_constants_no_user_defined_by_hand.cpp)
create_standard_test(NAME multiple_grid_cells SOURCES test_multiple_grid_cells.cpp)
if(ENABLE_JSON)
  create_standard_test(NAME rate_constants_no_user_defined_example_with_config SOURCES test_rate_constants_no_user_defined_with_config.cpp)
endif()

################################################################################
# Copy test data

add_custom_target(copy_tutorial_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/configs ${CMAKE_BINARY_DIR}/configs)

################################################################################
# Create zip folders of the configurations so that they can be 
# placed in the docs _static folder for downloads

# ~~~~~~~~~~~~~~~~
# There appears to be a bug with the archive create: https://gitlab.kitware.com/cmake/cmake/-/issues/25260
# ~~~~~~~~~~~~~~~~


# file(GLOB subdirectories LIST_DIRECTORIES true "${CMAKE_CURRENT_SOURCE_DIR}/configs/*")

# # Loop through the list of directories
# foreach(subdirectory ${subdirectories})
#     # Check if the current item is a directory
#     if(IS_DIRECTORY ${subdirectory})
#         # Extract the directory name from the path
#         get_filename_component(directory_name ${subdirectory} NAME)

#         # set the output file name
#         # set(output_archive "${CMAKE_SOURCE_DIR}/docs/source/_static/tutorials/${directory_name}.zip")
#         set(output_archive "${CMAKE_CURRENT_SOURCE_DIR}/configs/${directory_name}.zip")

#         message(STATUS "subdirectory: ${subdirectory}")

#         # create the archive
#         file(ARCHIVE_CREATE
#           OUTPUT "${output_archive}"
#           PATHS 
#               "${subdirectory}"
#           FORMAT "zip"
#         )

#         message(STATUS "Created archive: ${output_archive}")
#     endif()
# endforeach()