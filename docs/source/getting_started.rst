
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

.. code-block:: C++

  #include <iomanip>
  #include <iostream>
  #include <micm/process/arrhenius_rate_constant.hpp>
  #include <micm/solver/rosenbrock.hpp>

  using namespace micm;

  int main(const int argc, const char *argv[])
  {
    auto foo = Species{ "Foo" };
    auto bar = Species{ "Bar" };
    auto baz = Species{ "Baz" };

    Phase gas_phase{ std::vector<Species>{ foo, bar, baz } };

    System chemical_system{ SystemParameters{ .gas_phase_ = gas_phase } };

    Process r1 = Process::create()
                     .reactants({ foo })
                     .products({ Yield(bar, 0.8), Yield(baz, 0.2) })
                     .rate_constant(ArrheniusRateConstant({ .A_ = 1.0e-3 }))
                     .phase(gas_phase);

    Process r2 = Process::create()
                     .reactants({ foo, bar })
                     .products({ Yield(baz, 1) })
                     .rate_constant(ArrheniusRateConstant({ .A_ = 1.0e-5, .C_ = 110.0 }))
                     .phase(gas_phase);

    std::vector<Process> reactions{ r1, r2 };

    RosenbrockSolver solver{ chemical_system, reactions, RosenbrockSolverParameters{} };

    State state = solver.GetState();

    state.conditions_[0].temperature_ = 287.45;  // K
    state.conditions_[0].pressure_ = 101319.9;   // Pa
    state.SetConcentration(foo, 20.0);           // mol m-3

    std::cout << "foo,       bar,      baz" << std::endl;
    for (int i = 0; i < 10; ++i)
    {
      auto result = solver.Solve(500.0, state);
      state.variables_ = result.result_;
      std::cout << std::fixed << std::setprecision(6)
                << state.variables_[0][state.variable_map_["Foo"]] << ", "
                << state.variables_[0][state.variable_map_["Bar"]] << ", "
                << state.variables_[0][state.variable_map_["Baz"]] << std::endl;
    }

    return 0;
  }

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

