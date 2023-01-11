# Overwrite the init values choosen by CMake
if (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
   set(CMAKE_Fortran_FLAGS "-O2 -i4 -gopt -Mextend -byteswapio -Mflushz -Kieee -Mnofma")
endif()
