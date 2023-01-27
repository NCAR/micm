#!/bin/bash

# turn on command echoing
set -v
# move to the directory this script is in
cd ${0%/*}
# define a function for failure tests
failure_test () {
  local expected_failure=$(echo $1 | sed -n 's/\([[:digit:]]\+\).*/\1/p')
  local output=$(../../util_property_set $1 2>&1)
  local failure_code=$(echo $output | sed -n 's/[[:space:]]*ERROR (Musica-\([[:digit:]]\+\).*/\1/p')
  if ! [ "$failure_code" = "$expected_failure" ]; then
    echo "Expected failure $expected_failure"
    echo "Got output: $output"
    exit 1
  else
    echo $output
  fi
}

if ! ../../util_property_set; then
  exit 1
fi

failure_test 732626019
failure_test 934288374

exit 0
