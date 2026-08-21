window.BENCHMARK_DATA = {
  "lastUpdate": 1787283061305,
  "repoUrl": "https://github.com/NCAR/micm",
  "entries": {
    "TS1 Wall-Clock Timing": [
      {
        "commit": {
          "author": {
            "email": "kyle.shores44@gmail.com",
            "name": "Kyle Shores",
            "username": "K20shores"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "36d59591a047f0b78695977badf90b3b1beb7619",
          "message": "Ts1 benchmark (#1057)\n\n* Add TS1 benchmark\n\n- refactor the benchmark so that we can share it between TS1, chapman\n- generate the TS1 mechanism header with python from the config file\n- chart the TS1 benchmarks in gh-pages\n- Add MICM_ENABLE_BENCHMARK, a dependent option on MICM_ENABLE_TESTS that\n  defaults to OFF.\n- set GLYOXAL diffusion coeffient to the other defaults instead of 1e300\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-08-11T14:19:36-05:00",
          "tree_id": "564acc88fb6092878592ea1317296751dc377082",
          "url": "https://github.com/NCAR/micm/commit/36d59591a047f0b78695977badf90b3b1beb7619"
        },
        "date": 1786479140514,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 118282.77,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 117919.35,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 96999.57,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 82655.25,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 79572.83,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 47473.54,
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
          "id": "40f2019d70d50f51ab8132ef5f9af87eb94a725d",
          "message": "Auto-format code changes (#1058)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-11T15:35:40-05:00",
          "tree_id": "bc3e2a66d1e5a8d363a6b118b24992b9f856e873",
          "url": "https://github.com/NCAR/micm/commit/40f2019d70d50f51ab8132ef5f9af87eb94a725d"
        },
        "date": 1786483603950,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 114222.76,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 114775.63,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 91350.98,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 79425.12,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 75055.5,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 43065.93,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "kyle.shores44@gmail.com",
            "name": "Kyle Shores",
            "username": "K20shores"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "95fa70322b57bc4e22ac67099a35558a4e8fcf1a",
          "message": "Chart the CUDA solver on the GPU runner, not its CPU (#1059)\n\nChart the CUDA solver on the GPU runner\n\nThe CIRRUS a10 job called bench_chapman.sh and profile_chapman.sh with no\nbackend argument, and the scripts default to cpu. Explicitly pass gpu\nso that cirrus only reports gpu numbers\n\nAlso, remove valgrind for GPU runner because it cannot measure the\ninstruction count on the GPU and doesn't provide more information than\ncounting the cpu runs.\n\nCo-authored-by: Claude Opus 5 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-08-11T15:50:01-05:00",
          "tree_id": "03a0a93c80614b414d93501bdc476639d76d9acf",
          "url": "https://github.com/NCAR/micm/commit/95fa70322b57bc4e22ac67099a35558a4e8fcf1a"
        },
        "date": 1786484739570,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 123722.89,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 123160.2,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 101840.87,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 87949.66,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 86729.67,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 49781.02,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "kyle.shores44@gmail.com",
            "name": "Kyle Shores",
            "username": "K20shores"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6c1491ce289b51ef58d444a2b9176f43501868c3",
          "message": "adding more detail to performance table and description (#1060)",
          "timestamp": "2026-08-14T10:15:31-05:00",
          "tree_id": "e240f5661808519cbd26a09f6d8af87868f4e520",
          "url": "https://github.com/NCAR/micm/commit/6c1491ce289b51ef58d444a2b9176f43501868c3"
        },
        "date": 1786723675440,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 117830.11,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 116866.35,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 96752.69,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 83954.96,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 79892.12,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 46042.9,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "kyle.shores44@gmail.com",
            "name": "Kyle Shores",
            "username": "K20shores"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "97ac9e5d8aadd345c242722ee8274d71dfe0f73e",
          "message": "update python version to correct ci failure (#1061)",
          "timestamp": "2026-08-14T10:17:55-05:00",
          "tree_id": "7e9cd39e9e1beef97b5d3c7547b0b3dfd6cfe9b2",
          "url": "https://github.com/NCAR/micm/commit/97ac9e5d8aadd345c242722ee8274d71dfe0f73e"
        },
        "date": 1786723759027,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 115329.32,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 115160.81,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 92414.77,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 81495.91,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 77471.47,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 45142.36,
            "unit": "ms"
          }
        ]
      },
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
        "date": 1787273125519,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 85023.75,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 86526,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 83726.26,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 68690.29,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 66875.51,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 38603.83,
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
        "date": 1787283060585,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 113609.72,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 112977.2,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 99226.89,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 82565.19,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 74634.68,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 44447.7,
            "unit": "ms"
          }
        ]
      }
    ]
  }
}