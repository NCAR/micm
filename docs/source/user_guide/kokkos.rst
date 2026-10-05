.. _Kokkos:

Kokkos
======

MICM can build its solvers on top of `Kokkos <https://kokkos.org>`_, giving a single
source tree that runs on CPUs (Serial or OpenMP) and on GPUs (CUDA, HIP) through
Kokkos' performance-portability layer.  This is separate from MICM's hand-written
CUDA backend; the two are selected by different CMake options.

Enabling the backend
--------------------

Kokkos support is off by default.  Turn it on with ``MICM_ENABLE_KOKKOS`` and select a
Kokkos execution backend with Kokkos' own options:

.. code-block:: console

  # host, serial
  $ cmake -S . -B build -DMICM_ENABLE_KOKKOS=ON -DKokkos_ENABLE_SERIAL=ON

  # host, OpenMP
  $ cmake -S . -B build -DMICM_ENABLE_KOKKOS=ON -DKokkos_ENABLE_OPENMP=ON

  # NVIDIA GPU -- set the architecture for your device, and widen the vector length
  # (see "Choosing the vector width" below -- the default of 4 is far too narrow for a GPU)
  $ cmake -S . -B build -DMICM_ENABLE_KOKKOS=ON -DKokkos_ENABLE_CUDA=ON \
      -DKokkos_ARCH_AMPERE86=ON -DMICM_DEFAULT_VECTOR_SIZE=32

If ``MICM_ENABLE_KOKKOS=ON`` and no Kokkos backend is requested, Kokkos falls back to
Serial -- the build succeeds and runs on the host, so confirm the backend you asked for
is the one you got.  Architecture options are named ``Kokkos_ARCH_<ARCH>`` (for example
``Kokkos_ARCH_AMPERE86``, ``Kokkos_ARCH_HOPPER90``); an unrecognized name is silently
ignored by CMake and leaves the architecture to autodetection.

Kokkos is fetched at configure time if it is not already installed.  To build against
an existing installation, point CMake at it with ``-DCMAKE_PREFIX_PATH=/path/to/kokkos``
and it will be used instead of a fresh download.  This is the recommended route when you
intend to ``install`` MICM, because a Kokkos-enabled MICM package re-resolves Kokkos via
``find_dependency(Kokkos)`` and needs it to be discoverable.

.. note::
   Enabling Kokkos changes how *every* MICM header is compiled: it defines
   ``MICM_USE_KOKKOS``, which switches the ``MICM_LAMBDA``, ``MICM_DEVICE_FUNCTION`` and
   ``MICM_CONSTEXPR`` macros over to their Kokkos spellings.  Do not mix objects built
   with and without ``MICM_ENABLE_KOKKOS`` in one binary.

Linking
-------

.. code-block:: cmake

  find_package(micm REQUIRED)

  add_executable(my_target my_target.cpp)
  target_link_libraries(my_target PUBLIC musica::micm_kokkos)

``musica::micm_kokkos`` and ``musica::micm`` are equivalent in a Kokkos-enabled build;
both carry the Kokkos dependency and the ``MICM_USE_KOKKOS`` definition.

Initializing Kokkos
-------------------

Kokkos must be initialized before any MICM Kokkos type is constructed, and every such
object must be destroyed before Kokkos is finalized.  A ``Kokkos::View`` that outlives
``Kokkos::finalize()`` aborts the process.  The simplest correct form is a scope guard:

.. code-block:: cpp

  #include <micm/Kokkos.hpp>
  #include <Kokkos_Core.hpp>

  int main(int argc, char** argv)
  {
    Kokkos::ScopeGuard guard(argc, argv);

    // build solvers and states inside this scope
    return 0;
  }

Host and device synchronization
-------------------------------

The Kokkos matrices keep a host copy and a device copy, and MICM does **not**
synchronize them for you.  Writing to a state through its host accessors updates only
host memory; the solver reads device memory.  You must call ``CopyToDevice()`` after
writing inputs and ``CopyToHost()`` before reading results.

There is no dirty-flag tracking and no diagnostic if you forget, and on a Serial or
OpenMP build host and device memory are the same allocation -- so a missing copy is
invisible there and only shows up as wrong answers on a GPU build.

Worked example
--------------

.. code-block:: cpp

  #include <micm/Kokkos.hpp>

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  auto solver = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(gas_phase))
                    .SetReactions({ r1, r2, r3, r4 })
                    .Build();

  auto state = solver.GetState(n_grid_cells);

  // write inputs on the host ...
  state.variables_[0] = concentrations;
  state.custom_rate_parameters_[0] = photo_rates;
  state.conditions_[0].temperature_ = 284.19;
  state.conditions_[0].pressure_ = 101253.3;

  // ... then push them to the device
  state.variables_.CopyToDevice();
  state.custom_rate_parameters_.CopyToDevice();
  state.conditions_.CopyToDevice();

  solver.UpdateStateParameters(state);
  auto result = solver.Solve(30.0, state);

  // pull results back before reading them on the host
  state.variables_.CopyToHost();

Choosing the vector width
-------------------------

Like the CPU vectorized solver (:ref:`Vectorized matrix solver`), the Kokkos matrices are
blocked over ``L`` grid cells.  ``L`` is the second template parameter of
``micm::KokkosSolverBuilder`` and defaults to ``MICM_DEFAULT_VECTOR_SIZE`` (4).  That single
CMake cache variable sets the width for the CPU, CUDA and Kokkos matrix types alike -- there is
no separate Kokkos vector-size variable -- and ``micm::KokkosDenseMatrix<T>`` and
``micm::KokkosSparseMatrix<T>`` take the same default:

.. code-block:: cpp

  // 128 grid cells per block instead of the default 4
  using WideBuilder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters, 128>;

``L`` is also the width of the intra-team parallel loop, so on a GPU it sets how much of
each team has work to do: the team size is ``min(L, team_size_max)``.  A default of 4 therefore
means 4-thread CUDA blocks and leaves most of a warp idle, so a GPU build should raise it --
``-DMICM_DEFAULT_VECTOR_SIZE=32`` is one warp and a reasonable starting point.  Because the same
variable also sets the CPU width, note that ``docs/performance.md`` measures wider values as
faster on the host too, so raising it is not purely a GPU concern.
Note that per-thread temporaries in the solver scale with ``L``, so very large values
trade occupancy for local-memory pressure.  Benchmark your mechanism and grid size --
the best value is problem dependent.

Reproducibility
---------------

The Kokkos error norm is accumulated across teams with atomic floating-point additions,
so the summation order -- and therefore the adaptive step-size sequence -- is not
reproducible run to run, and does not match the CPU or CUDA backends bit for bit.  If
you need bit-reproducible results, use the CPU or CUDA backend.

Available types
---------------

``micm/Kokkos.hpp`` exports the Kokkos equivalents of the CPU aliases.  Each one is blocked at
``MICM_DEFAULT_VECTOR_SIZE``, the same width the CPU and CUDA aliases use:

.. list-table::
   :header-rows: 1

   * - Type
     - Description
   * - ``micm::KokkosSolverBuilder<Params, L>``
     - Builder for a Kokkos-backed solver
   * - ``micm::KokkosRosenbrock``
     - Rosenbrock solver at the default vector width
   * - ``micm::KokkosBackwardEuler``
     - Backward Euler solver at the default vector width
   * - ``micm::KokkosState``
     - State for the above
   * - ``micm::KokkosDenseReal`` / ``micm::KokkosSparseReal``
     - Underlying dense and sparse matrix types
