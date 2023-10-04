
Getting Started
===============

Build and Test
--------------

CPU
~~~
To build and install MICM locally, you must have the following libraries installed:

- [sphinx](https://github.com/sphinx-doc/sphinx)
- [sphinx-book-theme](https://github.com/executablebooks/sphinx-book-theme)
- [sphinx-design](https://github.com/executablebooks/sphinx-design)
- [breathe](https://github.com/breathe-doc/breathe)

You must also have CMake installed on your machine.

Open a terminal window, navigate to a folder where you would like the MICM files to exist,
and run the following commands::

    $ git clone https://github.com/NCAR/micm.git
    $ cd micm
    $ mkdir build
    $ cd build
    $ ccmake ..
    $ make install -j 8
    $ make test

CMake will allow for setting options such as the installation directory
with CMAKE_INSTALL_PREFIX, or various build flags such as BUILD_DOCS, ENABLE_CUDA, etc.

Docker Container
~~~~~~~~~~~~~~~~

Build and run the image::

    $ docker build -t micm -f Dockerfile.nvhpc .
    $ docker run --rm -it micm

If you would like, you can ssh into a running docker container and edit the files there.

GPU
~~~

NCAR Hardware
-------------

On Cheyenne
~~~~~~~~~~~

On Casper
~~~~~~~~~

On Gust and Derecho
~~~~~~~~~~~~~~~~~~~
To compile and test on gust::

    $ qinteractive -A NTDD0005 --ngpus=1
    $ module load cmake/3.25.2 nvhpc/23.1 cuda/11.7.1
    $ mkdir build && cd build
    $ cmake -DENABLE_OPENACC=OFF -DENABLE_CUDA=ON -D GPU_TYPE="a100" ..
    $ make
    $ make test

NOAA Hardware
-------------

Run an Example
--------------
The following example solves the fictitious chemical system::

  foo       --k1--> 0.8 bar + 0.2 baz
  foo + bar --k2--> baz

The `k1` and `k2` rate constants are for Arrhenius reactions.
See the [MICM documentation](https://ncar.github.io/micm/)
for details on the types of reactions available in MICM and how to configure them.
To solve this system save the following code in a file named `foo_chem.cpp`

.. literalinclude:: ../../test/tutorial/test_README_example.cpp
  :language: cpp

To build and run the example using GNU::

  $ g++ -o foo_chem foo_chem.cpp -I<CMAKE_INSTALL_PREFIX>/include -std=c++20
  $ ./foo_chem

You should see an output including this

.. csv-table::
   :header: "time [s]", "foo", "bar", "baz"
   :widths: 6, 10, 10, 10

    0.000000,  11.843503,   5.904845,   1.907012
    500.000000,   6.792023,   9.045965,   3.317336
    1000.000000,   3.828700,  10.740589,   4.210461
    1500.000000,   2.138145,  11.663685,   4.739393
    2000.000000,   1.187934,  12.169452,   5.042503
    2500.000000,   0.658129,  12.447502,   5.213261
    3000.000000,   0.363985,  12.600676,   5.308597
    3500.000000,   0.201076,  12.685147,   5.361559
    4000.000000,   0.111028,  12.731727,   5.390884
    4500.000000,   0.061290,  12.757422,   5.407096

