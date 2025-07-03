.. _JIT:

A JIT-compiled solver
=====================

Solving a chemical mechanism is a computationally intensive task. There are many reasons for this such as
multiple iterations by a solver to achieve an integration, stiff systems requiring internal substepping to 
acheive a numerically stable solution, and cache misses, among others. 

This tutorial focuses on alleviating the cache misses. A popular method for handling cache misses is to pre-compute the indices. 
This method, which may be referred to as ahead-of-time (AOT) compilation, is used in applications such as KPP :cite:`Damian2002`.
Pre-computed methods require code preprocessors and preclude runtime configurable software, which is a goal of MICM.

MICM uses just-in-time (JIT) compiled functions built with `LLVM JIT <https://llvm.org/docs/tutorial/BuildingAJIT1.html>`_ libraries
to supply runtime-configurable chemistry to avoid cache misses in important chemistry functions.

Up until now, a :cpp:class:`micm::RosenbrockSolver` has been used. This is a special class in MICM which builds all
of the componenets needed to solve chemistry in memory, including the forcing function and the jacobian. MICM
also provides a :cpp:class:`micm::JitRosenbrockSolver` which builds and compiles several important chemistry functions at runtime.

What are we JIT-compilng?
-------------------------

A list of compiled functions is below. How they are compiled and the ellided operations are beyond the scope of this tutorial, but you are free to inspect the source code yourself.

- :cpp:func:`micm::RosenbrockSolver::AlphaMinusJacobian` is JIT compiled and called with :cpp:func:`micm::JitRosenbrockSolver::AlphaMinusJacobian`
- :cpp:func:`micm::LuDecomposition::Decompose` is JIT compiled and called with :cpp:func:`micm::JitLuDecomposition::Decompose`
- :cpp:func:`micm::LinearSolver::Factor` and :cpp:func:`micm::LinearSolver::Solve` are JIT compiled and called with :cpp:func:`micm::JitLinearSolver::Factor` and :cpp:func:`micm::JitLinearSolver::Solve`
- :cpp:func:`micm::ProcessSet::AddForcingTerms` and :cpp:func:`micm::ProcessSet::AddJacobianTerms` are JIT compiled and called with :cpp:func:`micm::JitProcessSet::AddForcingTerms` and :cpp:func:`micm::JitProcessSet::AddJacobianTerms`

So what does this gain me?
--------------------------

Runtime configuraiton of chemical mechanisms that are *fast*. Let's compare the speed of :cpp:class:`micm::RosenbrockSolver` 
and the :cpp:class:`micm::JitRosenbrockSolver` classes applied to the same problem. We will use the simple
fictitous chemical system from previous examples for simplicity, but feel free to try out more complex
mechanisms to see the effects of JIT compiling on compute time in more realistic cases.

If you're looking for a copy and paste, choose
the appropriate tab below and be on your way! Otherwise, stick around for a line by line explanation.


.. tab-set::

    .. tab-item:: JIT-compiled Rosenbrock Solver

        .. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
          :language: cpp

To build and run the example using GNU (assuming the default install location), copy and
paste the example code into a file named ``foo_jit_chem.cpp`` and run:

.. code:: bash

  g++ -o foo_jit_chem foo_jit_chem.cpp -I/usr/local/micm-3.2.0/include `llvm-config --cxxflags --ldflags --system-libs --libs support core orcjit native irreader` -std=c++20 -fexceptions
  ./foo_jit_chem

Line-by-line explanation
------------------------

Starting with the header files, we need headers for timing, output, and of course for both types of solvers.

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 1-5

Next, we use our namespace, define our number of gridcells (1 for now), and some
partial template specializations. We are using our custom vectorized matrix, which
groups mutliple reactions across grid cells into tiny blocks in a vector, allowing multiple
grid cells to be solved simultaneously.

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 7-16

Now, all at once, is the function which runs either type of solver. We set all species
concentrations to 1 :math:`\mathrm{mol\ m^-3}`.
Additionally, we are collecting all of the solver stats across all solving timesteps

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 18-71

Finally, the main function which reads the configuration and initializes the jit solver.

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 73-113

The only additional step here is to make an instance of the :cpp:class:`micm::JitCompiler`
and pass it as a shared pointer to the :cpp:class:`micm::JitRosenbrockSolver`. We also are
using our vectorized matrix for both solvers. The :cpp:class:`micm::JitRosenbrockSolver` only works
with the vectorized matrix whereas the :cpp:class:`micm::RosenbrockSolver` works with a regular matrix.
At construction of the :cpp:class:`micm::JitRosenbrockSolver`, all JIT functions are compiled. We record that
time here.

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 103-113

Finally, we run both solvers, output their cumulative stats, and compare their results.


.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 115-181


The output will be similar to this:

.. code:: bash

  Jit compile time: 38305334 nanoseconds
  Standard solve time: 122582 nanoseconds
  JIT solve time: 96541 nanoseconds
  Standard solve stats: 
          accepted: 130
          function_calls: 260
          jacobian_updates: 130
          number_of_steps: 130
          accepted: 130
          rejected: 0
          decompositions: 130
          solves: 390
          singular: 0
          total_update_state_time: 416 nanoseconds
          total_forcing_time: 14163 nanoseconds
          total_jacobian_time: 9965 nanoseconds
          total_linear_factor_time: 24893 nanoseconds
          total_linear_solve_time: 12842 nanoseconds

  JIT solve stats: 
          accepted: 130
          function_calls: 260
          jacobian_updates: 130
          number_of_steps: 130
          accepted: 130
          rejected: 0
          decompositions: 130
          solves: 390
          singular: 0
          total_update_state_time: 458 nanoseconds
          total_forcing_time: 5036 nanoseconds
          total_jacobian_time: 2168 nanoseconds
          total_linear_factor_time: 17958 nanoseconds
          total_linear_solve_time: 7586 nanoseconds