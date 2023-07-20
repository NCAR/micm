MICM Chemistry
==============

Model Independent Chemical Mechanisms.

[![License](https://img.shields.io/github/license/NCAR/micm.svg)](https://github.com/NCAR/micm/blob/master/LICENSE)
[![CI Status](https://github.com/NCAR/micm/actions/workflows/test.yml/badge.svg)](https://github.com/NCAR/micm/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/NCAR/micm/branch/main/graph/badge.svg?token=ATGO4DKTMY)](https://codecov.io/gh/NCAR/micm)

Copyright (C) 2018-2020 National Center for Atmospheric Research


## Build and run

### CPU

### GPU

#### Checking that it compiles on your local machine

1. Build the image
```
docker build -t micm -f Dockerfile.nvhpc .
```
2. Run the container
```
docker run --rm -it micm
```
3. Compile micm. After running the previous command, you can run `make` and see your compile errors.
```
make
```
4. If you'd like, you can ssh into a running docker container and edit the files there.

#### On Gust
```
qinteractive -A NTDD0005 --ngpus=1
module load cmake/3.25.2 nvhpc/23.1 cuda/11.7.1
mkdir build && cd build
cmake -DENABLE_OPENACC=OFF -DENABLE_CUDA=ON -D GPU_TYPE="a100" ..
make
make test
```

# Building and installing
To build and install MICM locally, you must have the following libraries installed:

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
docker run -it ghcr.io/ncar/micm:release bash
```

Inside the container, you can run the MICM tests from the `/build/` folder:

```
cd build/
# to run the tests
make test
```

## Local installation

```
git clone https://github.com/NCAR/micm.git
cd micm
mkdir build
cd build
ccmake ..
make install -j 8
# to run the tests
make test
```

## Install

After completing the previous step run `sudo make install`.
This will install the MICM static library.

If you would later like to uninstall MICM, you can run
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

- [Apache 2.0](/LICENSE)

Copyright (C) 2018-2023 National Center for Atmospheric Research
