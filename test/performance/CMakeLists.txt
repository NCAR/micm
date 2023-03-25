################################################################################
# Link micmlib to a test and add it to the suite

macro(add_std_test test_name)
  target_include_directories(${test_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_name} micmlib)
  if (ENABLE_NSYS)
     add_test(NAME ${test_name} COMMAND nsys profile --force-overwrite true -o nsys_out --trace openacc,cuda ${CMAKE_BINARY_DIR}/${test_name})
  else()
     if (ENABLE_MPI)
        if (ENABLE_OPENACC)
            add_test(NAME ${test_name} COMMAND mpirun -v -np ${NUM_TASKS} ${CMAKE_BINARY_DIR}/set_device_rank.sh ${CMAKE_BINARY_DIR}/${test_name})
        else()
            add_test(NAME ${test_name} COMMAND mpirun -v -np ${NUM_TASKS} ${CMAKE_BINARY_DIR}/${test_name})
        endif()
     else()
        add_test(NAME ${test_name} COMMAND ${test_name})
     endif()
  endif()
endmacro(add_std_test)

################################################################################
# MICM performance tests

add_executable(test_performance driver.F90)
if (ENABLE_OPENACC)
   if (CMAKE_Fortran_COMPILER STREQUAL "NVFORTRAN")
      set_source_files_properties(driver.F90 PROPERTIES COMPILE_FLAGS "-acc -gpu=cc70,lineinfo,nofma -Minfo=all")
      set_target_properties(test_performance PROPERTIES LINK_FLAGS "-acc -gpu=cc70,lineinfo,nofma -Minfo=all")
   else()    # GNU
      set_source_files_properties(driver.F90 PROPERTIES COMPILE_FLAGS "-fopenacc")
      set_target_properties(test_performance PROPERTIES LINK_FLAGS "-fopenacc")
   endif()
endif()
add_std_test(test_performance)

################################################################################
