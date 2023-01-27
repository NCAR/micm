
musica-core
===========

A library of model components and utilities.

[![License](https://img.shields.io/github/license/NCAR/musica-core)](https://github.com/NCAR/musica-core/blob/main/LICENSE)
[![CI Status](https://github.com/NCAR/musica-core/actions/workflows/test.yml/badge.svg)](https://github.com/NCAR/musica-core/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/NCAR/musica-core/branch/main/graph/badge.svg?token=WIBA0JE3OE)](https://codecov.io/gh/NCAR/musica-core)

Copyright (C) 2020 National Center for Atmospheric Research

A working draft of the `musica-core` documentation can be found [here](https://ncar.github.io/musica-core/html/index.html).

# Install the library and run tests

The `musica-core` library requires:

- git
- a Fortran compiler
- CMake ([https://cmake.org/download/](https://cmake.org/download/))
- NetCDF ([https://www.unidata.ucar.edu/downloads/netcdf/](https://www.unidata.ucar.edu/downloads/netcdf/))
- json-fortran ([https://github.com/jacobwilliams/json-fortran](https://github.com/jacobwilliams/json-fortran))

To install the library, run (you may need to set cmake options if the dependencies are not automatically found):

```
git clone --recurse-submodules https://github.com/NCAR/musica-core
cd musica-core
mkdir build
cd build
cmake ..
make
```

To run the library test, from the `build/` folder run:

```
make test
```
