# Overwrite the init values choosen by CMake
if (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
   set(CMAKE_Fortran_FLAGS "-O2")
   set(CMAKE_EXE_LINKER_FLAGS "-O2 -acc -gpu=cc70,lineinfo,nofma -Minfo=all")
endif()
