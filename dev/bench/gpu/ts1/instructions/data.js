window.BENCHMARK_DATA = {
  "lastUpdate": 1786479810210,
  "repoUrl": "https://github.com/NCAR/micm",
  "entries": {
    "TS1 Instruction Counts (GPU runner, CPU backend)": [
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
        "date": 1786479809434,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 42204327615,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 42261153385,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 28030365028,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 21945909115,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 18436654904,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 14653787145,
            "unit": "instructions"
          }
        ]
      }
    ]
  }
}