window.BENCHMARK_DATA = {
  "lastUpdate": 1786723676187,
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
      }
    ]
  }
}