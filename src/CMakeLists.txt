################################################################################
# micm source

set(MICM_SRC micm.cpp)

add_library(micmlib ${MICM_SRC})
target_include_directories(micmlib PUBLIC ${CMAKE_SOURCE_DIR}/include)

################################################################################
