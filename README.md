MICM Chemistry
==============

Model Independent Chemical Mechanisms.

[![License](https://img.shields.io/github/license/NCAR/micm.svg)](https://github.com/NCAR/micm/blob/master/LICENSE)
[![CI Status](https://github.com/NCAR/micm/actions/workflows/test.yml/badge.svg)](https://github.com/NCAR/micm/actions/workflows/test.yml)


# Building and installing
To build and install MICM locally, you must have the following libraries installed:

- [json-fortran](https://github.com/jacobwilliams/json-fortran)

- [sphinx](https://github.com/sphinx-doc/sphinx)
- [sphinx-book-theme](https://github.com/executablebooks/sphinx-book-theme)
- [sphinx-design](https://github.com/executablebooks/sphinx-design)
- [breathe](https://github.com/breathe-doc/breathe)

You must also have CMake installed on your machine. 

To install MICM locally,
open a terminal window, navigate to a folder where you would like the MICM files to exist,
and run the following commands:

## Build and run (Docker version)

To build and run the stand-alone version of MICM, you must have [Docker Desktop](https://www.docker.com/get-started) installed and running. With Docker Desktop running, open a terminal window and run the following command to start the MICM container:

```
docker run -it ghcr.io/ncar/mcim:release bash
```

Inside the container, you can run the MICM tests from the `/build/` folder:

```
cd build/
# to run the tests
make test
# to use the standalone tool
./micm examples/full_config.json
```

## Build and run (local build version)

```
git clone https://github.com/NCAR/micm.git
cd mcim
mkdir build
cd build
ccmake ..
make -j 8
```

You will now have a runnable exectubable for `mcim` and the tests in the build directory.

`./mcim examples/full_config.json`.

Inspect the output file `photolysis_rate_constants.nc` to see the results!

## Install

After completing the previous step run `sudo make install`.
This wll install the tuvx static library, the MCIM configuration
and runtime data, as well as the standalone `mcim` exectuable, which can be
added to your system path to make the executable useable from any directory.

If you would later lake to uninstall tuv-x, you can run
`sudo make uninstall` from the build directory.

# Citation


# Community and contributions
We welcome contributions and feedback from anyone, everything from updating
the content or appearance of the documentation to new and
cutting edge science.

- [Collaboration](https://github.com/NCAR/musica/blob/main/docs/Software%20Development%20Plan.pdf)
  - Anyone interested in scientific collaboration
which would add new software functionality should read the [MUSICA software development plan](https://github.com/NCAR/musica/blob/main/docs/Software%20Development%20Plan.pdf).

- [Code of conduct]
  - Please read this through to you understand the expectations with how to interact with this project.

- [Contributor's guide]
  - Before submiitting a PR, please thouroughly read this to you understand our expectations. We reserve the right to reject any PR not meeting our guidelines.


# Documentation
Please see the [MICM documentation](https://ncar.github.io/micm/) for detailed
installation and usage instructions.

# License

- [GPL 2.0](/LICENSE)

Copyright (C) 2018-2020 National Center for Atmospheric Research
