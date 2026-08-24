window.BENCHMARK_DATA = {
  "lastUpdate": 1787608073842,
  "repoUrl": "https://github.com/NCAR/micm",
  "entries": {
    "Chapman Timing (Kokkos backend)": [
      {
        "commit": {
          "author": {
            "email": "bbakernoaa@users.noreply.github.com",
            "name": "Barry Baker",
            "username": "bbakernoaa"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "359d392c96e6acd1b02c2d1af8cbf4b4c6a19f75",
          "message": "Implement Kokkos support for MICM chemical kinetics and optimizations (#967)\n\n* feat: Add Kokkos support for chemical forcing and Jacobian terms\n\nThis commit implements a Kokkos backend for MICM, following the design of\nthe existing CUDA driver. It includes:\n- KokkosDenseMatrix and KokkosSparseMatrix using Kokkos Views\n- KokkosProcessSet with vectorized kernels for forcing and Jacobian\n- Integration into the CMake build system (MICM_ENABLE_KOKKOS)\n- Unit tests for the new Kokkos components\n\n* feat(kokkos): Add SetAlgebraicVariableIds method and optimize forcing/Jacobian indexing\n\n- Add SetAlgebraicVariableIds method to KokkosProcessSet for setting algebraic variable flags on device\n- Optimize forcing term indexing to use number_of_forcing_species instead of number_of_species\n- Optimize Jacobian term indexing by pre-computing jacobian_group_size and jac_group_offset to reduce redundant calculations\n- Refactor Jacobian index computation to use pre-computed offset for better performance\n- Update test file with copyright header and consolidate test utilities using shared test policy\n- Improve memory access patterns in AddForcingTerms and SubtractJacobianTerms kernels\n\n* style(kokkos): Format code for improved readability\n\n- Break long lines in loop assignments to improve code readability\n- Split long Kokkos::View declarations across multiple lines\n- Reformat assignment statements in SetJacobianFlatIds method\n- Improve line length consistency throughout kokkos_process_set.hpp\n- Update corresponding test file formatting to match header changes\n\n* ci(workflows): Add Kokkos CUDA GPU testing job\n\n- Add new kokkos_cuda_tests job to runner workflow\n- Configure NVIDIA GPU container with CUDA 13.0 support\n- Install required dependencies (cmake, g++)\n- Build with Kokkos CUDA enabled on Ampere architecture\n- Run full test suite on GPU hardware with parallel execution\n- Enable GPU resource allocation via Cirrus group configuration\n\n* chore(dependencies): Update Kokkos dependency from 4.4.01 to 4.7.04\n\n- Bump Kokkos version to 4.7.04 for improved performance and bug fixes\n- Ensures compatibility with latest Kokkos features and optimizations\n\n* refactor(kokkos): Replace auto with explicit template parameters in process set\n\n- Replace `auto` with explicit `StatePolicy` template parameter in `AddForcingTerms` method\n- Replace `auto` with explicit `StatePolicy` template parameter in `SubtractJacobianTerms` method\n- Improve code clarity and type safety by making state parameter types explicit\n- Enhance IDE support and type checking for Kokkos process operations\n\n* chore(kokkos): Update to Kokkos 5.1.0 and improve code formatting\n\n- Update Kokkos dependency from 4.7.04 to 5.1.0 in cmake/dependencies.cmake\n- Replace HostMirror type alias with explicit Kokkos::View<T*, Kokkos::HostSpace> in dense and sparse matrix headers\n- Remove RandomAccess memory trait from ProcessSetParam view declarations for improved compatibility\n- Add blank lines between include groups for better code organization\n- Reorder includes alphabetically (type_traits before vector)\n- Add braces to single-statement for loops for consistency and readability\n- Format GroupVectorSize() methods across multiple lines for improved readability\n\n* ci(workflows): Add Kokkos serial backend testing and refactor utilities\n\n- Add Kokkos serial backend CI jobs to macOS, Ubuntu, and Windows workflows\n- Remove kokkos_util.hpp header and consolidate utility includes in main Kokkos.hpp\n- Refactor KokkosProcessSet to remove redundant Initialize() calls and simplify device struct updates\n- Clean up trailing whitespace in workflow files\n- Reorganize include order in Kokkos.hpp for better dependency management\n- Update test CMakeLists.txt to remove kokkos_util.hpp references\n- Improve code formatting and consistency across Kokkos utility headers\n\n* refactor(kokkos): Optimize matrix synchronization with unmanaged host views\n\n- Remove persistent h_view_ member from KokkosDenseMatrix and KokkosSparseMatrix\n- Replace manual element-by-element loops with Kokkos::deep_copy for better performance\n- Use unmanaged host views in CopyToDevice() and CopyToHost() to avoid redundant allocations\n- Simplify constructor initialization by removing h_view_ creation\n- Add comprehensive documentation explaining device/host synchronization pattern\n- Improve Fill() method by removing unnecessary h_view_ management\n- Reduces memory overhead and aligns with CUDA matrix implementation pattern\n\n* test(kokkos): Add comprehensive unit tests for KokkosSparseMatrix\n\n- Add test_kokkos_sparse_matrix.cpp with full test coverage for KokkosSparseMatrix\n- Test zero matrix initialization and const zero matrix behavior\n- Test scalar assignment and diagonal operations\n- Test single and multi-block matrix operations with various ordering policies\n- Test print and print non-zero functionality across all vector orderings\n- Add using declaration for base class assignment operator in KokkosSparseMatrix\n- Register new test executable in CMakeLists.txt with proper linking and test configuration\n- Ensures KokkosSparseMatrix implementation is fully validated across different matrix configurations\n\n* test(kokkos): Expand KokkosSparseMatrix unit tests with BlockFunction coverage\n\n- Rename type aliases from KokkosVectorOrdering to KokkosOrdering for clarity\n- Add BlockFunction infrastructure tests (ArrayFunction, MultiMatrixArrayFunction)\n- Add temporary variable and view reuse tests (MultipleTemporaries, BlockViewReuse)\n- Add function reusability and multi-matrix structure tests\n- Add sparse + dense/vector matrix interaction tests\n- Add ordering compatibility and error condition tests\n- Organize tests into logical sections with comments for maintainability\n- Expand test coverage across all four vector orderings (1-4) where applicable\n\n* remove process set; add dockerfile; debug tests\n\n* remove process info structs\n\n* add kokkos view classes and no-op funcs for matrices\n\n* add kokkos dense matrix functions\n\n* add kokkos sparse matrix functions\n\n* add reducer functions\n\n* full vector support in matrix functions\n\n* add missing template parameter\n\n* create serial and cuda Kokkos docker files\n\n* update for cuda lambda reqs\n\n* update lambda syntax for kokkos\n\n* update dense matrix Function lambdas and add Kokkos tests\n\n* add copy and fill tests for Kokkos dense matrix\n\n* add reduce sum kokkos dense matrix support\n\n* add kokkos dense matrix support for remaining reducers\n\n* add remaining Kokkos dense matrix tests\n\n* add dense matrix vector capture test\n\n* updates to sparse matrices\n\n* add template arg to LU decomp classes\n\n* debug actions and cuda code\n\n* debug action failures\n\n* exclude unsupported cuda tests\n\n* debug runner action failure\n\n* add LU mozart in place Kokkos support and tests\n\n* add kokkos support and tests to all LU decomps\n\n* add kokkos support and tests to linear solvers\n\n* add Kokkos support and tests for process set\n\n* add kokkos support for rate constants\n\n* add kokkos support and tests for constraints\n\n* add kokkos support and tests for all solvers\n\n* move scalar views to state and debug a couple test failures\n\n* move index structs to types header and fix cuda rate constant function builds\n\n* debug cuda build; add kokkos serial docker action\n\n* debug action failures\n\n* debug action failures\n\n* debug actions\n\n* add clang tidy suppressions and debug benchmark action\n\n* clang format\n\n* tidy\n\n* tidier\n\n* tidiest\n\n* debug tidy action\n\n* most tidiest\n\n* exclude umbrella headers in clang-tidy scans\n\n* continue to debug clang tidy failure\n\n* actually tidy, for real this time\n\n* use lambda for solver clamping function\n\n* move host view to class member for kokkos matrices\n\n* add Kokkos solvers to benchmark and create Kokkos runner action\n\n* update comments\n\n* bump nvhpc image version in kokkos runner and try to sidestep compiler error\n\n* reduce parallel build cores for kokkos runner\n\n* address review comments\n\n* remove unneeded vector allocation and double division\n\n* remove unneeded structs\n\n* remove *Const*  view classes\n\n* remove unneeded vector and scalar views\n\n* address review comments\n\n* remove commented out tests\n\n* address review comments\n\n* address Claude suggestions\n\n* address review comments\n\n* add Kokkos::AUTO back to team policy constructors\n\n---------\n\nCo-authored-by: bbakernoaa <22104759+bbakernoaa@users.noreply.github.com>\nCo-authored-by: Matt Dawson <matt@cohere-llc.com>",
          "timestamp": "2026-08-20T17:06:40-07:00",
          "tree_id": "93c3dcc5a2964761dfc7a8e18b3ae7fc57954e3a",
          "url": "https://github.com/NCAR/micm/commit/359d392c96e6acd1b02c2d1af8cbf4b4c6a19f75"
        },
        "date": 1787273907565,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 4272.13,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 2293.79,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 1307.6,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 755.53,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 581.94,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f94aeddd622bba69cffe955449c0c60c2691ba7c",
          "message": "Auto-format code changes (#1071)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-20T20:39:44-06:00",
          "tree_id": "75f5bb50d306515e6114d38bf5a320d834b3c1ce",
          "url": "https://github.com/NCAR/micm/commit/f94aeddd622bba69cffe955449c0c60c2691ba7c"
        },
        "date": 1787282959478,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 737.47,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 427.35,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 282.14,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 216.67,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 245.45,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "336ae6400e5b1f1810f987f3b0595a86263fc073",
          "message": "Auto-format code changes (#1077)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-21T13:16:08-06:00",
          "tree_id": "be7fa93762ea5f85a8a025e85e5c2a32f0190290",
          "url": "https://github.com/NCAR/micm/commit/336ae6400e5b1f1810f987f3b0595a86263fc073"
        },
        "date": 1787344482930,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 724.61,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 1574.6,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 1285.45,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 1154.48,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 1188.71,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "sunjian@ucar.edu",
            "name": "Jian Sun",
            "username": "sjsprecious"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e2e21ed02e38f42a5368752a28f9ed32affc9b83",
          "message": "Fix build, CI and packaging for Kokkos deps (#1075)\n\n* fix build, ci and packaging\n\n* fix the broken CI tests",
          "timestamp": "2026-08-21T14:55:43-06:00",
          "tree_id": "87ed9a7dc6ed9621c2458bc09e1421cebf3a6bd9",
          "url": "https://github.com/NCAR/micm/commit/e2e21ed02e38f42a5368752a28f9ed32affc9b83"
        },
        "date": 1787348919902,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 731.01,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 415.72,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 265.7,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 199.01,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 224.65,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "92f00b12a364fa38aa06578a4e3a591ba37ccaa0",
          "message": "Auto-format code changes (#1079)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-24T13:42:40-06:00",
          "tree_id": "c59d49e5a010597c9540d86e8b9802cfda565163",
          "url": "https://github.com/NCAR/micm/commit/92f00b12a364fa38aa06578a4e3a591ba37ccaa0"
        },
        "date": 1787608072865,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 722.72,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 409.03,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 264.17,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 197.99,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 223.65,
            "unit": "ms"
          }
        ]
      }
    ],
    "TS1 Timing (Kokkos backend)": [
      {
        "commit": {
          "author": {
            "email": "bbakernoaa@users.noreply.github.com",
            "name": "Barry Baker",
            "username": "bbakernoaa"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "359d392c96e6acd1b02c2d1af8cbf4b4c6a19f75",
          "message": "Implement Kokkos support for MICM chemical kinetics and optimizations (#967)\n\n* feat: Add Kokkos support for chemical forcing and Jacobian terms\n\nThis commit implements a Kokkos backend for MICM, following the design of\nthe existing CUDA driver. It includes:\n- KokkosDenseMatrix and KokkosSparseMatrix using Kokkos Views\n- KokkosProcessSet with vectorized kernels for forcing and Jacobian\n- Integration into the CMake build system (MICM_ENABLE_KOKKOS)\n- Unit tests for the new Kokkos components\n\n* feat(kokkos): Add SetAlgebraicVariableIds method and optimize forcing/Jacobian indexing\n\n- Add SetAlgebraicVariableIds method to KokkosProcessSet for setting algebraic variable flags on device\n- Optimize forcing term indexing to use number_of_forcing_species instead of number_of_species\n- Optimize Jacobian term indexing by pre-computing jacobian_group_size and jac_group_offset to reduce redundant calculations\n- Refactor Jacobian index computation to use pre-computed offset for better performance\n- Update test file with copyright header and consolidate test utilities using shared test policy\n- Improve memory access patterns in AddForcingTerms and SubtractJacobianTerms kernels\n\n* style(kokkos): Format code for improved readability\n\n- Break long lines in loop assignments to improve code readability\n- Split long Kokkos::View declarations across multiple lines\n- Reformat assignment statements in SetJacobianFlatIds method\n- Improve line length consistency throughout kokkos_process_set.hpp\n- Update corresponding test file formatting to match header changes\n\n* ci(workflows): Add Kokkos CUDA GPU testing job\n\n- Add new kokkos_cuda_tests job to runner workflow\n- Configure NVIDIA GPU container with CUDA 13.0 support\n- Install required dependencies (cmake, g++)\n- Build with Kokkos CUDA enabled on Ampere architecture\n- Run full test suite on GPU hardware with parallel execution\n- Enable GPU resource allocation via Cirrus group configuration\n\n* chore(dependencies): Update Kokkos dependency from 4.4.01 to 4.7.04\n\n- Bump Kokkos version to 4.7.04 for improved performance and bug fixes\n- Ensures compatibility with latest Kokkos features and optimizations\n\n* refactor(kokkos): Replace auto with explicit template parameters in process set\n\n- Replace `auto` with explicit `StatePolicy` template parameter in `AddForcingTerms` method\n- Replace `auto` with explicit `StatePolicy` template parameter in `SubtractJacobianTerms` method\n- Improve code clarity and type safety by making state parameter types explicit\n- Enhance IDE support and type checking for Kokkos process operations\n\n* chore(kokkos): Update to Kokkos 5.1.0 and improve code formatting\n\n- Update Kokkos dependency from 4.7.04 to 5.1.0 in cmake/dependencies.cmake\n- Replace HostMirror type alias with explicit Kokkos::View<T*, Kokkos::HostSpace> in dense and sparse matrix headers\n- Remove RandomAccess memory trait from ProcessSetParam view declarations for improved compatibility\n- Add blank lines between include groups for better code organization\n- Reorder includes alphabetically (type_traits before vector)\n- Add braces to single-statement for loops for consistency and readability\n- Format GroupVectorSize() methods across multiple lines for improved readability\n\n* ci(workflows): Add Kokkos serial backend testing and refactor utilities\n\n- Add Kokkos serial backend CI jobs to macOS, Ubuntu, and Windows workflows\n- Remove kokkos_util.hpp header and consolidate utility includes in main Kokkos.hpp\n- Refactor KokkosProcessSet to remove redundant Initialize() calls and simplify device struct updates\n- Clean up trailing whitespace in workflow files\n- Reorganize include order in Kokkos.hpp for better dependency management\n- Update test CMakeLists.txt to remove kokkos_util.hpp references\n- Improve code formatting and consistency across Kokkos utility headers\n\n* refactor(kokkos): Optimize matrix synchronization with unmanaged host views\n\n- Remove persistent h_view_ member from KokkosDenseMatrix and KokkosSparseMatrix\n- Replace manual element-by-element loops with Kokkos::deep_copy for better performance\n- Use unmanaged host views in CopyToDevice() and CopyToHost() to avoid redundant allocations\n- Simplify constructor initialization by removing h_view_ creation\n- Add comprehensive documentation explaining device/host synchronization pattern\n- Improve Fill() method by removing unnecessary h_view_ management\n- Reduces memory overhead and aligns with CUDA matrix implementation pattern\n\n* test(kokkos): Add comprehensive unit tests for KokkosSparseMatrix\n\n- Add test_kokkos_sparse_matrix.cpp with full test coverage for KokkosSparseMatrix\n- Test zero matrix initialization and const zero matrix behavior\n- Test scalar assignment and diagonal operations\n- Test single and multi-block matrix operations with various ordering policies\n- Test print and print non-zero functionality across all vector orderings\n- Add using declaration for base class assignment operator in KokkosSparseMatrix\n- Register new test executable in CMakeLists.txt with proper linking and test configuration\n- Ensures KokkosSparseMatrix implementation is fully validated across different matrix configurations\n\n* test(kokkos): Expand KokkosSparseMatrix unit tests with BlockFunction coverage\n\n- Rename type aliases from KokkosVectorOrdering to KokkosOrdering for clarity\n- Add BlockFunction infrastructure tests (ArrayFunction, MultiMatrixArrayFunction)\n- Add temporary variable and view reuse tests (MultipleTemporaries, BlockViewReuse)\n- Add function reusability and multi-matrix structure tests\n- Add sparse + dense/vector matrix interaction tests\n- Add ordering compatibility and error condition tests\n- Organize tests into logical sections with comments for maintainability\n- Expand test coverage across all four vector orderings (1-4) where applicable\n\n* remove process set; add dockerfile; debug tests\n\n* remove process info structs\n\n* add kokkos view classes and no-op funcs for matrices\n\n* add kokkos dense matrix functions\n\n* add kokkos sparse matrix functions\n\n* add reducer functions\n\n* full vector support in matrix functions\n\n* add missing template parameter\n\n* create serial and cuda Kokkos docker files\n\n* update for cuda lambda reqs\n\n* update lambda syntax for kokkos\n\n* update dense matrix Function lambdas and add Kokkos tests\n\n* add copy and fill tests for Kokkos dense matrix\n\n* add reduce sum kokkos dense matrix support\n\n* add kokkos dense matrix support for remaining reducers\n\n* add remaining Kokkos dense matrix tests\n\n* add dense matrix vector capture test\n\n* updates to sparse matrices\n\n* add template arg to LU decomp classes\n\n* debug actions and cuda code\n\n* debug action failures\n\n* exclude unsupported cuda tests\n\n* debug runner action failure\n\n* add LU mozart in place Kokkos support and tests\n\n* add kokkos support and tests to all LU decomps\n\n* add kokkos support and tests to linear solvers\n\n* add Kokkos support and tests for process set\n\n* add kokkos support for rate constants\n\n* add kokkos support and tests for constraints\n\n* add kokkos support and tests for all solvers\n\n* move scalar views to state and debug a couple test failures\n\n* move index structs to types header and fix cuda rate constant function builds\n\n* debug cuda build; add kokkos serial docker action\n\n* debug action failures\n\n* debug action failures\n\n* debug actions\n\n* add clang tidy suppressions and debug benchmark action\n\n* clang format\n\n* tidy\n\n* tidier\n\n* tidiest\n\n* debug tidy action\n\n* most tidiest\n\n* exclude umbrella headers in clang-tidy scans\n\n* continue to debug clang tidy failure\n\n* actually tidy, for real this time\n\n* use lambda for solver clamping function\n\n* move host view to class member for kokkos matrices\n\n* add Kokkos solvers to benchmark and create Kokkos runner action\n\n* update comments\n\n* bump nvhpc image version in kokkos runner and try to sidestep compiler error\n\n* reduce parallel build cores for kokkos runner\n\n* address review comments\n\n* remove unneeded vector allocation and double division\n\n* remove unneeded structs\n\n* remove *Const*  view classes\n\n* remove unneeded vector and scalar views\n\n* address review comments\n\n* remove commented out tests\n\n* address review comments\n\n* address Claude suggestions\n\n* address review comments\n\n* add Kokkos::AUTO back to team policy constructors\n\n---------\n\nCo-authored-by: bbakernoaa <22104759+bbakernoaa@users.noreply.github.com>\nCo-authored-by: Matt Dawson <matt@cohere-llc.com>",
          "timestamp": "2026-08-20T17:06:40-07:00",
          "tree_id": "93c3dcc5a2964761dfc7a8e18b3ae7fc57954e3a",
          "url": "https://github.com/NCAR/micm/commit/359d392c96e6acd1b02c2d1af8cbf4b4c6a19f75"
        },
        "date": 1787276205834,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 291247.66,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 162658.05,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 94450.5,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 52783.5,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 33183.46,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f94aeddd622bba69cffe955449c0c60c2691ba7c",
          "message": "Auto-format code changes (#1071)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-20T20:39:44-06:00",
          "tree_id": "75f5bb50d306515e6114d38bf5a320d834b3c1ce",
          "url": "https://github.com/NCAR/micm/commit/f94aeddd622bba69cffe955449c0c60c2691ba7c"
        },
        "date": 1787283414306,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 52056.88,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 29399.02,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 17451.33,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 11595.39,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 12856.09,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "336ae6400e5b1f1810f987f3b0595a86263fc073",
          "message": "Auto-format code changes (#1077)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-21T13:16:08-06:00",
          "tree_id": "be7fa93762ea5f85a8a025e85e5c2a32f0190290",
          "url": "https://github.com/NCAR/micm/commit/336ae6400e5b1f1810f987f3b0595a86263fc073"
        },
        "date": 1787345513899,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 120979.9,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 67522.08,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 39324.68,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 25655.78,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 25740.69,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "sunjian@ucar.edu",
            "name": "Jian Sun",
            "username": "sjsprecious"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e2e21ed02e38f42a5368752a28f9ed32affc9b83",
          "message": "Fix build, CI and packaging for Kokkos deps (#1075)\n\n* fix build, ci and packaging\n\n* fix the broken CI tests",
          "timestamp": "2026-08-21T14:55:43-06:00",
          "tree_id": "87ed9a7dc6ed9621c2458bc09e1421cebf3a6bd9",
          "url": "https://github.com/NCAR/micm/commit/e2e21ed02e38f42a5368752a28f9ed32affc9b83"
        },
        "date": 1787349377068,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 52076.13,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 29519.8,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 17608.03,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 11660.77,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 12814.8,
            "unit": "ms"
          }
        ]
      }
    ]
  }
}