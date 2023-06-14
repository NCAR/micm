MICM Chemistry
==============

Model Independent Chemical Mechanisms.

[![License](https://img.shields.io/github/license/NCAR/micm.svg)](https://github.com/NCAR/micm/blob/master/LICENSE)
[![CI Status](https://github.com/NCAR/micm/actions/workflows/test.yml/badge.svg)](https://github.com/NCAR/micm/actions/workflows/test.yml)

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
qinteractive -A<ACCOUNT_NUMBER> --ngpus=1
module load cmake/3.25.2 nvhpc/23.1 cuda/11.7.1
mkdir build && cd build
cmake -DENABLE_OPENACC=ON -DENABLE_CUDA=ON -DENABLE_REGESSION_TESTS=OFF ..
make
make test
```