window.BENCHMARK_DATA = {
  "lastUpdate": 1786483602469,
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
      }
    ]
  }
}