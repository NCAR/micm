window.BENCHMARK_DATA = {
  "lastUpdate": 1786131713533,
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
          "message": "Add tool for identifying possible regressions (#1053)\n\n* draft performance regression tests\n\n* update tests; us continuous benchmarking action\n\n* debug cuda benchmark\n\n* Potential fix for pull request finding\n\nCo-authored-by: Copilot Autofix powered by AI <175728472+Copilot@users.noreply.github.com>\n\n* Potential fix for pull request finding\n\nCo-authored-by: Copilot Autofix powered by AI <175728472+Copilot@users.noreply.github.com>\n\n* Potential fix for pull request finding\n\nCo-authored-by: Copilot Autofix powered by AI <175728472+Copilot@users.noreply.github.com>\n\n* Potential fix for pull request finding\n\nCo-authored-by: Copilot Autofix powered by AI <175728472+Copilot@users.noreply.github.com>\n\n* Potential fix for pull request finding\n\nCo-authored-by: Copilot Autofix powered by AI <175728472+Copilot@users.noreply.github.com>\n\n* address copilot comments\n\n* update actions\n\n* update actions\n\n* update actions\n\n* debug actions\n\n* debug actions\n\n* Install valgrind on the GPU runner and stop hiding script failures\n\nThe nvhpc container image does not ship valgrind, so\nscripts/profile_chapman.sh exited at its guard and printed nothing. The\nprofile steps used a pipeline, and the container runs steps under `sh -e`,\nwhich has no pipefail. The step reported success, and\nto_benchmark_json.py wrote an empty JSON array. The benchmark action then\nfailed with a misleading message about a missing benchmark result.\n\n- Install valgrind before the build, so chapman_bench also finds\n  <valgrind/callgrind.h> and scopes the instruction count to Solve().\n- Write each script's output to a file, then convert that file. A plain\n  command failure stops the step under `sh -e`, with no shell dependency.\n- Upload the raw TSV output, to make a future failure easy to diagnose.\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>\n\n* Give the GPU runner its own benchmark series and fix the data paths\n\nThe comment steps read dev/bench, which is the action default, but the store\nsteps wrote dev/bench/instructions and dev/bench/timing. The timing series\nalso carried two names: \"Chapman Timing\" in the comment step and \"Chapman\nWall-Clock Timing\" in the store step. github-action-benchmark keys its\nhistory by name and directory, so no comparison could ever match. The step\nstill reports success, which hides the problem.\n\nrunner.yml and benchmark-charts.yml also wrote the same two series names to\nthe same two directories on a push to main, both with auto-push. The two\nmachines differ: runner.yml measures the CIRRUS a10 GPU runner, and\nbenchmark-charts.yml measures ubuntu-latest. One shared series mixes them.\n\n- Move the GPU runner data to dev/bench/gpu/{instructions,timing}.\n- Tag its series names with \"(GPU runner)\".\n- Use the same name and directory in each comment step and its store step.\n\nThe disjoint directories also settle the concurrent push to gh-pages, since\nthe action retries a rejected push with a rebase and the two workflows now\ntouch different files.\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>\n\n* Point the documentation links at Read the Docs and link the benchmark charts\n\nThe docs build on Read the Docs now, and RTD builds every commit. The gh-pages\nsite stopped updating on 2026-02-07, but five links still sent readers there.\n\n- Repoint the five links to https://micm.readthedocs.io/en/latest/.\n- Add the trailing underscore in getting_started.rst. Without it the reST\n  target is interpreted text, so that link never worked.\n- Add a Performance section to the README with the four chart links. The table\n  names the machine for each chart, because the GPU runner and ubuntu-latest\n  keep separate series and must not be compared against each other.\n\nThe chart URLs return 404 until the first push to main creates the data.\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>\n\n* Point the intersphinx mappings at Read the Docs\n\nAll three sibling projects moved to Read the Docs. The MechanismConfiguration\ntarget was already dead, because that repository does not enable GitHub Pages.\nEach new inventory returns 200.\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>\n\n---------\n\nCo-authored-by: Copilot Autofix powered by AI <175728472+Copilot@users.noreply.github.com>\nCo-authored-by: Kyle Shores <kyle.shores44@gmail.com>\nCo-authored-by: Claude Opus 5 (1M context) <noreply@anthropic.com>",
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
      }
    ]
  }
}