################################################################################
# Tests

create_standard_test(NAME jit_linear_solver SOURCES test_jit_linear_solver.cpp)
create_standard_test(NAME jit_lu_decomposition_doolittle SOURCES test_jit_lu_decomposition_doolittle.cpp)
create_standard_test(NAME jit_rosenbrock SOURCES test_jit_rosenbrock.cpp)
create_standard_test(NAME jit_solver_builder SOURCES test_jit_solver_builder.cpp)