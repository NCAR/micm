window.BENCHMARK_DATA = {
  "lastUpdate": 1786723757608,
  "repoUrl": "https://github.com/NCAR/micm",
  "entries": {
    "TS1 Instruction Counts": [
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
        "date": 1786479138273,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 42091902505,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 42297342364,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 28138463861,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 22005221525,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 18514826722,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 14645751568,
            "unit": "instructions"
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
        "date": 1786483601623,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 42091902505,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 42297342364,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 28138463861,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 22005221525,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 18514826722,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 14645751568,
            "unit": "instructions"
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
        "date": 1786484737944,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 42091902505,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 42297342364,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 28138463861,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 22005221525,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 18514826722,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 14645751568,
            "unit": "instructions"
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
        "date": 1786723673119,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 42091902505,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 42297342364,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 28138463861,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 22005221525,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 18514826722,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 14645751568,
            "unit": "instructions"
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
        "date": 1786723756723,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 42091902505,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 42297342364,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 28138463861,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 22005221525,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 18514826722,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 14645751568,
            "unit": "instructions"
          }
        ]
      }
    ]
  }
}