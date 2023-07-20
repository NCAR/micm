MICM Chemistry
==============

Model Independent Chemical Module.

[![License](https://img.shields.io/github/license/NCAR/micm.svg)](https://github.com/NCAR/micm/blob/master/LICENSE)
[![CI Status](https://github.com/NCAR/micm/actions/workflows/test.yml/badge.svg)](https://github.com/NCAR/micm/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/NCAR/micm/branch/main/graph/badge.svg?token=ATGO4DKTMY)](https://codecov.io/gh/NCAR/micm)

Copyright (C) 2018-2020 National Center for Atmospheric Research

# Getting Started

## Installing MICM locally
To build and install MICM locally, you must have the following libraries installed:

- [sphinx](https://github.com/sphinx-doc/sphinx)
- [sphinx-book-theme](https://github.com/executablebooks/sphinx-book-theme)
- [sphinx-design](https://github.com/executablebooks/sphinx-design)
- [breathe](https://github.com/breathe-doc/breathe)

You must also have CMake installed on your machine. 

Open a terminal window, navigate to a folder where you would like the MICM files to exist,
and run the following commands:

```
git clone https://github.com/NCAR/micm.git
cd micm
mkdir build
cd build
ccmake ..
sudo make install -j 8
```

To run the tests:

```
make test
```

If you would later like to uninstall MICM, you can run
`sudo make uninstall` from the `build/` directory.

## Running a MICM Docker container

You must have [Docker Desktop](https://www.docker.com/get-started) installed and running. With Docker Desktop running, open a terminal window and run the following command to start the MICM container:

```
docker run -it ghcr.io/ncar/micm:release bash
```

Inside the container, you can run the MICM tests from the `/build/` folder:

```
cd /build/
make test
```


# Citation

MICM is part of the MUSICA project and can be cited by reference to the MUSICA vision paper. The BibTeX entry below can be used to generate a citation for this.

```
@Article { acom.software.musica-vision,
    author = "Gabriele G. Pfister and Sebastian D. Eastham and Avelino F. Arellano and Bernard Aumont and Kelley C. Barsanti and Mary C. Barth and Andrew Conley and Nicholas A. Davis and Louisa K. Emmons and Jerome D. Fast and Arlene M. Fiore and Benjamin Gaubert and Steve Goldhaber and Claire Granier and Georg A. Grell and Marc Guevara and Daven K. Henze and Alma Hodzic and Xiaohong Liu and Daniel R. Marsh and John J. Orlando and John M. C. Plane and Lorenzo M. Polvani and Karen H. Rosenlof and Allison L. Steiner and Daniel J. Jacob and Guy P. Brasseur",
    title = "The Multi-Scale Infrastructure for Chemistry and Aerosols (MUSICA)",
    journal = "Bulletin of the American Meteorological Society",
    year = "2020",
    publisher = "American Meteorological Society",
    address = "Boston MA, USA",
    volume = "101",
    number = "10",
    doi = "10.1175/BAMS-D-19-0331.1",
    pages= "E1743 - E1760",
    url = "https://journals.ametsoc.org/view/journals/bams/101/10/bamsD190331.xml"
}
```

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
