window.BENCHMARK_DATA = {
  "lastUpdate": 1786131715416,
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
      }
    ]
  }
}