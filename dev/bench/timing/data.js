window.BENCHMARK_DATA = {
  "lastUpdate": 1786720775707,
  "repoUrl": "https://github.com/NCAR/micm",
  "entries": {
    "Chapman Wall-Clock Timing": [
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
        "date": 1786129124251,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1230.24,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1239.38,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 812.93,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 574.98,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 508.91,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 420.4,
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
        "date": 1786131714982,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1092.29,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1094.77,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 752.42,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 589.18,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 523.87,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 471.24,
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
        "date": 1786476231156,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1254.13,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1268.86,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 796.86,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 556.72,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 492.74,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 409.86,
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
        "date": 1786480789266,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1260.41,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1268.98,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 781.1,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 559.97,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 492.85,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 407.14,
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
        "date": 1786481661811,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1216.91,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1226.6,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 802.36,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 582.03,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 523.05,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 431.75,
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
        "date": 1786720774967,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 1263.07,
            "unit": "ms"
          },
          {
            "name": "vector1",
            "value": 1272.67,
            "unit": "ms"
          },
          {
            "name": "vector2",
            "value": 780.71,
            "unit": "ms"
          },
          {
            "name": "vector4",
            "value": 562.54,
            "unit": "ms"
          },
          {
            "name": "vector8",
            "value": 497.04,
            "unit": "ms"
          },
          {
            "name": "vector128",
            "value": 413.17,
            "unit": "ms"
          }
        ]
      }
    ]
  }
}