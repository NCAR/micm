######################################################################
# MICM library

set(MICM_SRC core.F90)

add_library(micmlib ${MICM_SRC})

target_link_libraries(micmlib musicalib)

set_target_properties(micmlib PROPERTIES OUTPUT_NAME micm)
