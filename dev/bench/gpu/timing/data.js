window.BENCHMARK_DATA = {
  "lastUpdate": 1788191693615,
  "repoUrl": "https://github.com/NCAR/micm",
  "entries": {
    "Chapman Timing (GPU runner)": [
      {
        "commit": {
          "author": {
            "email": "mattldawson@gmail.com",
            "name": "Matt Dawson",
            "username": "mattldawson"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "38411382debc5fd7419fc0d294f8630963c1bff6",
          "message": "Add tool for identifying possible regressions (#1053)",
          "timestamp": "2026-08-07T13:55:29-05:00",
          "tree_id": "a5025974f3fdc6f15268c8379d6cf7987aa5c3ed",
          "url": "https://github.com/NCAR/micm/commit/38411382debc5fd7419fc0d294f8630963c1bff6"
        },
        "date": 1786130202486,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1159.58,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1158.01,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 806.63,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 611.8,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 549.25,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 487.66,
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
          "id": "0e3f1d1299cf456228d88988c6492f5642355c55",
          "message": "Auto-format code changes (#1056)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-07T12:38:51-07:00",
          "tree_id": "0808a0402701e22a6b98506c0488a75572e0c6c4",
          "url": "https://github.com/NCAR/micm/commit/0e3f1d1299cf456228d88988c6492f5642355c55"
        },
        "date": 1786132335176,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1215.43,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1165.07,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 787.86,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 627.18,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 563.74,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 497.42,
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
          "id": "36d59591a047f0b78695977badf90b3b1beb7619",
          "message": "Ts1 benchmark (#1057)\n\n* Add TS1 benchmark\n\n- refactor the benchmark so that we can share it between TS1, chapman\n- generate the TS1 mechanism header with python from the config file\n- chart the TS1 benchmarks in gh-pages\n- Add MICM_ENABLE_BENCHMARK, a dependent option on MICM_ENABLE_TESTS that\n  defaults to OFF.\n- set GLYOXAL diffusion coeffient to the other defaults instead of 1e300\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-08-11T14:19:36-05:00",
          "tree_id": "564acc88fb6092878592ea1317296751dc377082",
          "url": "https://github.com/NCAR/micm/commit/36d59591a047f0b78695977badf90b3b1beb7619"
        },
        "date": 1786476872311,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1165.63,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1172.03,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 791.59,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 610.42,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 567.44,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 496.11,
            "unit": "ms"
          }
        ]
      }
    ],
    "Chapman Timing (GPU backend)": [
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
        "date": 1786482230969,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 114.93,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 74.62,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 61.59,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 61.41,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 56.35,
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
        "date": 1786721482207,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 116.58,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 73.51,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 62.14,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 60.81,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 56.34,
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
        "date": 1787271627078,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 214.74,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 152.37,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 119.75,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 111.62,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 106.89,
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
        "date": 1787281295735,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 114.66,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 78.37,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 68.19,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 68.01,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 67.71,
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
        "date": 1787343348072,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 113.11,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 86.57,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 81.29,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 81.23,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 81,
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
        "date": 1787346883391,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 113.78,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 77.91,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 72.46,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 72.68,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 71.86,
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
        "date": 1787601558947,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 115.43,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 78.78,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 67.3,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 67.57,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 66.86,
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
          "id": "0b489c8a6098665e03a5e129a3c2aedd22b06e1e",
          "message": "Auto-format code changes (#1086)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-27T09:12:35-06:00",
          "tree_id": "8e37e9c287aac11bfadfce2b73e4c96b712d4d57",
          "url": "https://github.com/NCAR/micm/commit/0b489c8a6098665e03a5e129a3c2aedd22b06e1e"
        },
        "date": 1787844497963,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 214.31,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 152.2,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 119.48,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 111.79,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 106.94,
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
          "id": "7c62ea7691ddbd816b06ee7622727d3a31840d66",
          "message": "Fix the broken GPU benchmark PR comments (#1087)\n\nfix GPU pr comments",
          "timestamp": "2026-08-28T10:45:48-06:00",
          "tree_id": "daad4e7006aa865820e9778cf635a86b7ff8b681",
          "url": "https://github.com/NCAR/micm/commit/7c62ea7691ddbd816b06ee7622727d3a31840d66"
        },
        "date": 1787936358135,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 214.91,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 152.55,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 119.77,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 111.78,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 105.35,
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
          "id": "bb291a5f284c23a5579e43d36ae4683bac27c614",
          "message": "Fix the occasionally broken single-precision CI test (#1088)\n\nfix sp ci failure",
          "timestamp": "2026-08-31T09:28:43-06:00",
          "tree_id": "72dcbf43a2ace27e93cdae40f5cf3407ff6290de",
          "url": "https://github.com/NCAR/micm/commit/bb291a5f284c23a5579e43d36ae4683bac27c614"
        },
        "date": 1788191023020,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 212.93,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 151.9,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 119.81,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 111.77,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 106.95,
            "unit": "ms"
          }
        ]
      }
    ],
    "TS1 Timing (GPU backend)": [
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
        "date": 1786482416354,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 21956.75,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 10979.25,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 6941.37,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 5024.98,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 4215.15,
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
        "date": 1786721667612,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 22001.05,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 11004.34,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 6964.77,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 5059.25,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 4260.8,
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
        "date": 1787272416231,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 111351.51,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 51072.68,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 24960.26,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 13158.99,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 9118.88,
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
        "date": 1787281514687,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 26195.46,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 13178.78,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 8305.05,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 5959.25,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 4990.87,
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
        "date": 1787343609237,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 31481.07,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 15757.27,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 9913.85,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 7106.35,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 5954.17,
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
        "date": 1787347121431,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 28339.42,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 14230.91,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 8946.86,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 6398.23,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 5385.82,
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
        "date": 1787601779536,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 26372.32,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 13192.06,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 8342.15,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 6010.1,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 5042.53,
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
          "id": "0b489c8a6098665e03a5e129a3c2aedd22b06e1e",
          "message": "Auto-format code changes (#1086)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-27T09:12:35-06:00",
          "tree_id": "8e37e9c287aac11bfadfce2b73e4c96b712d4d57",
          "url": "https://github.com/NCAR/micm/commit/0b489c8a6098665e03a5e129a3c2aedd22b06e1e"
        },
        "date": 1787845290923,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 120745.57,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 49950.15,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 24419.61,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 13254.51,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 9191.2,
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
          "id": "7c62ea7691ddbd816b06ee7622727d3a31840d66",
          "message": "Fix the broken GPU benchmark PR comments (#1087)\n\nfix GPU pr comments",
          "timestamp": "2026-08-28T10:45:48-06:00",
          "tree_id": "daad4e7006aa865820e9778cf635a86b7ff8b681",
          "url": "https://github.com/NCAR/micm/commit/7c62ea7691ddbd816b06ee7622727d3a31840d66"
        },
        "date": 1787937144538,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 112704.98,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 49965.25,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 24397.69,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 13261.92,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 9198.67,
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
          "id": "bb291a5f284c23a5579e43d36ae4683bac27c614",
          "message": "Fix the occasionally broken single-precision CI test (#1088)\n\nfix sp ci failure",
          "timestamp": "2026-08-31T09:28:43-06:00",
          "tree_id": "72dcbf43a2ace27e93cdae40f5cf3407ff6290de",
          "url": "https://github.com/NCAR/micm/commit/bb291a5f284c23a5579e43d36ae4683bac27c614"
        },
        "date": 1788191692870,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "vector1",
            "value": 102697.18,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 42427.57,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 20687.66,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 11237.15,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 7770.12,
            "unit": "ms"
          }
        ]
      }
    ]
  }
}