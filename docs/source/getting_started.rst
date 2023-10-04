###############
Getting Started
###############

Build and Test
==============

CPU
---
To build and install MICM locally, you must have the following libraries installed:

- `CMake <https://cmake.org/>`_
  - `installation <https://cmake.org/download/>`_

Then, it's enough for you to configure and install micm on your computer. Because micm is header-only library, the install
step will simply copy the header files into the normal location required by your system.

.. code-block:: console
  
    $ git clone https://github.com/NCAR/micm.git
    $ cd micm
    $ mkdir build
    $ cd build
    $ ccmake ..
    $ make install -j 8
    $ make test

CMake will allow for setting options such as the installation directory
with CMAKE_INSTALL_PREFIX, or various build flags such as BUILD_DOCS, ENABLE_CUDA, etc.

MICM can optionally include support for json configuration reading, OpenMP,
JIT-compiled chemistry functions, and GPUs. Each of these requires an additional library. 
Some of these libraries can be included automatically with cmake build options,
others require that you have libraries installed on your system.

- JSON configuration support
  - When building micm, you need to enable the JSON option. This will download and configure the `nlohmann/jsoncpp library <https://github.com/nlohmann/json>`_  for you. For example: ``cmake -DENABLE_JSON=ON ..``
- JIT-compiled chemistry functions 
  - This requires `LLVM <https://llvm.org/docs/index.html>`_ to be installed with on your system. Once it is, you can include the jit options with ``cmake -DENBABLE_LLVM=ON ..``
- GPU support
  - Coming soon
- OpenMP
  - On macOS, you either need to configure cmake to use gcc which ships with OpenMP (either ``CXX=g++ cmake -DENABLE_OPENMP=ON ..`` or ``cmake -DCMAKE_CXX_COMPILER=g++ -DENABLE_OPENMP=ON ..``)

Docker Container
----------------

Build and run the image::

    $ docker build -t micm -f Dockerfile.nvhpc .
    $ docker run --rm -it micm

If you would like, you can ssh into a running docker container and edit the files there.

GPU
---

NCAR Hardware
^^^^^^^^^^^^^

On Cheyenne
^^^^^^^^^^^

On Casper
^^^^^^^^^

On Gust and Derecho
^^^^^^^^^^^^^^^^^^^

To compile and test on gust::

    $ qinteractive -A NTDD0005 --ngpus=1
    $ module load cmake/3.25.2 nvhpc/23.1 cuda/11.7.1
    $ mkdir build && cd build
    $ cmake -DENABLE_OPENACC=OFF -DENABLE_CUDA=ON -D GPU_TYPE="a100" ..
    $ make
    $ make test

NOAA Hardware
^^^^^^^^^^^^^

Run an Example
--------------
The following example solves the fictitious chemical system::

  foo       --k1--> 0.8 bar + 0.2 baz
  foo + bar --k2--> baz

The `k1` and `k2` rate constants are for Arrhenius reactions.
See the `MICM documentation <https://ncar.github.io/micm/>`
for details on the types of reactions available in MICM and how to configure them.
To solve this system save the following code in a file named `foo_chem.cpp`

.. literalinclude:: ../../test/tutorial/test_README_example.cpp
  :language: cpp

To build and run the example using GNU::

  $ g++ -o foo_chem foo_chem.cpp -I<CMAKE_INSTALL_PREFIX>/include -std=c++20
  $ ./foo_chem

Output::

  foo,       bar,      baz
  19.034389, 0.762719, 0.197464
  18.105748, 1.478520, 0.395242
  17.213802, 2.150391, 0.592159
  16.358113, 2.781130, 0.787212
  15.538111, 3.373351, 0.979560
  14.753106, 3.929496, 1.168498
  14.002317, 4.451851, 1.353446
  13.284884, 4.942548, 1.533932
  12.599887, 5.403583, 1.709582
  11.946359, 5.836817, 1.880104

