
musica-core
===========

A library of model components and utilities.

[![License](https://img.shields.io/github/license/NCAR/musica-core)](https://github.com/NCAR/musica-core/blob/master/LICENSE) [![Build Status](https://travis-ci.com/NCAR/musica-core.svg?branch=master)](https://travis-ci.com/NCAR/musica-core)

Copyright (C) 2020 National Center for Atmospheric Research

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

# Documentation

The musica-core documentation can be built using [Doxygen](https://www.doxygen.nl). After [installing](https://www.doxygen.nl/download.html) Doxygen, from the root `musica-core/` folder run:

```
cd doc
doxygen
```
Then, open `musica-core/doc/html/index.html` in a browser.

The documentation includes detailed instructions for using the `musica-core` library to build a model.

