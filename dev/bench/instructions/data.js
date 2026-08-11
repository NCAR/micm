window.BENCHMARK_DATA = {
  "lastUpdate": 1786480787741,
  "repoUrl": "https://github.com/NCAR/micm",
  "entries": {
    "Chapman Instruction Counts": [
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
        "date": 1786129122522,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 538524426,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 539264413,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 362547603,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 287543280,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 243199638,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 211176086,
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
          "id": "0e3f1d1299cf456228d88988c6492f5642355c55",
          "message": "Auto-format code changes (#1056)\n\nAuto-format code using Clang-Format\n\nCo-authored-by: GitHub Actions <actions@github.com>",
          "timestamp": "2026-08-07T12:38:51-07:00",
          "tree_id": "0808a0402701e22a6b98506c0488a75572e0c6c4",
          "url": "https://github.com/NCAR/micm/commit/0e3f1d1299cf456228d88988c6492f5642355c55"
        },
        "date": 1786131713038,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 538524426,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 539264413,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 362547603,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 287543280,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 243199638,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 211176086,
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
          "id": "36d59591a047f0b78695977badf90b3b1beb7619",
          "message": "Ts1 benchmark (#1057)\n\n* Add TS1 benchmark\n\n- refactor the benchmark so that we can share it between TS1, chapman\n- generate the TS1 mechanism header with python from the config file\n- chart the TS1 benchmarks in gh-pages\n- Add MICM_ENABLE_BENCHMARK, a dependent option on MICM_ENABLE_TESTS that\n  defaults to OFF.\n- set GLYOXAL diffusion coeffient to the other defaults instead of 1e300\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-08-11T14:19:36-05:00",
          "tree_id": "564acc88fb6092878592ea1317296751dc377082",
          "url": "https://github.com/NCAR/micm/commit/36d59591a047f0b78695977badf90b3b1beb7619"
        },
        "date": 1786476228688,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 530969733,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 539574531,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 366258762,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 287643944,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 247372360,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 211095607,
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
        "date": 1786480787014,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "standard",
            "value": 530969733,
            "unit": "instructions"
          },
          {
            "name": "vector1",
            "value": 539574531,
            "unit": "instructions"
          },
          {
            "name": "vector2",
            "value": 366258762,
            "unit": "instructions"
          },
          {
            "name": "vector4",
            "value": 287643944,
            "unit": "instructions"
          },
          {
            "name": "vector8",
            "value": 247372360,
            "unit": "instructions"
          },
          {
            "name": "vector128",
            "value": 211095607,
            "unit": "instructions"
          }
        ]
      }
    ]
  }
}