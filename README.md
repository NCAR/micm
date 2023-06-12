MICM Chemistry
==============

Model Independent Chemical Mechanisms.

[![License](https://img.shields.io/github/license/NCAR/micm.svg)](https://github.com/NCAR/micm/blob/master/LICENSE)
[![CI Status](https://github.com/NCAR/micm/actions/workflows/test.yml/badge.svg)](https://github.com/NCAR/micm/actions/workflows/test.yml)

Copyright (C) 2018-2020 National Center for Atmospheric Research


## Build and run

### CPU

### GPU

```
mkdir build && cd build
cmake -DENABLE_GPU=ON -DENABLE_REGESSION_TESTS=OFF ..
make
make test
```