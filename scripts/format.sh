#!/usr/bin/env bash
# Reformat the same files the Clang-Format workflow targets. Keep in sync
# with .github/workflows/clang_format_non_inline.yml.

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root"

CLANG_FORMAT="${CLANG_FORMAT:-clang-format}"

find include -type f \( -name '*.hpp' -o -name '*.h' -o -name '*.cuh' -o -name '*.inl' \) \
  -exec "$CLANG_FORMAT" -i --style=file --verbose {} +
find src -type f \( -name '*.cu' -o -name '*.hpp' -o -name '*.h' -o -name '*.cpp' \) \
  -exec "$CLANG_FORMAT" -i --style=file --verbose {} +
find test -type f \( -name '*.hpp' -o -name '*.h' -o -name '*.cpp' -o -name '*.cuh' -o -name '*.cu' \) \
  ! -path 'test/tutorial/*' \
  -exec "$CLANG_FORMAT" -i --style=file --verbose {} +
