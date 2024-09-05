# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-src")
  file(MAKE_DIRECTORY "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-src")
endif()
file(MAKE_DIRECTORY
  "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-build"
  "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix"
  "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/tmp"
  "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src/yaml-cpp-populate-stamp"
  "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src"
  "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src/yaml-cpp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src/yaml-cpp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/mthind/Desktop/code/micm/build/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src/yaml-cpp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
