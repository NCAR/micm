.. _JIT:

A JIT-compiled solver
=====================

Solving a chemical mechanism is a computationally intensive task. There are many reasons for this such as
multiple iterations by a solver to achieve an integration, stiff systems requiring internal substepping to 
acheive a numerically stable solution, and cache misses, among others. 

This tutorial focuses on alleviating the cache misses. A popular method for handling cache misses is to pre-compute the indices. 
This method, which may be referred to as ahead-of-time (AOT) compilation, is used in applications such as KPP :cite:`Damian2002`.
Pre-computed methods require code preprocessors and preclude runtime configurable software, which is a goal of micm.

MICM uses just-in-time (JIT) compiled functions built with `LLVM JIT <https://llvm.org/docs/tutorial/BuildingAJIT1.html>`_ libraries
to supply runtime-configurable chemistry to avoid cache misses in important chemistry functions.

Up until now, a :cpp:class:`micm::RosenbrockSolver` has been used. This is a special class in micm which builds all
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
and the :cpp:class:`micm::JitRosenbrockSolver` classes applied to the same problem. This time we'll use the carbon-bond 5
mechanism :cite:`Yarwood2005`, which is an update to the carbon bond mechanism from :cite:`Gery1989`.

If you're looking for a copy and paste, choose
the appropriate tab below and be on your way! Otherwise, stick around for a line by line explanation.


.. tab-set::

    .. tab-item:: OpenAtmos Configuration reading

        .. raw:: html

            <div class="download-div">
              <a href="../_static/tutorials/carbon_bond_5.zip" download>
                <button class="download-button">Download zip configuration</button>
              </a>
            </div>

        .. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
          :language: cpp

Line-by-line explanation
------------------------

Starting with the header files, we need headers for timing, outpu, reading the configuration,
and of course for both types of solvers.

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 1-7

Next, we use our namespace, define our number of gridcells (1 for now), and some
partial template specializations. We are using our custom vectorized matrix, which
groups mutliple reactions across grid cells into tiny blocks in a vector, allowing multiple
grid cells to be solved simultaneously.

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 14-18

Now, all at once, is the function which runs either type of solver. We set all species
concentrations to 1 :math:`\mathrm{mol m^-3}` and set the rate paramters for all of the
photolysis and emissions reactions. Additionally, we are collecting all of the solver
stats across all solving timesteps

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 20-96

Finally, the main function which reads the configuration and initializes the jit solver.

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 101-117

The only additionall step here is to make an instance of the :cpp:class:`micm::JitCompiler`
and pass it as a shared pointer to the :cpp:class:`micm::JitRosenbrockSolver`. We also are
using our vectorized matrix for both solvers. The :cpp:class:`micm::JitRosenbrockSolver` only works
with the vectorized matrix whereas the :cpp:class:`micm::RosenbrockSolver` works with a regular matrix.
At construction of the :cpp:class:`micm::JitRosenbrockSolver`, all JIT functions are compiled. We record that
time here.

.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 119-140

Finally, we run both solvers, output their cumulative stats, and compare their results.


.. literalinclude:: ../../../test/tutorial/test_jit_tutorial.cpp
  :language: cpp
  :lines: 142-198


The output will be similar to this:

.. code:: bash

  Jit compile time: 207422864917 nanoseconds
  Result time: 8593750 nanoseconds
  JIT result time: 4904666 nanoseconds
  Result stats: 
  accepted: 162
  function_calls: 334
  jacobian_updates: 162
  number_of_steps: 172
  accepted: 162
  rejected: 0
  decompositions: 172
  solves: 516
  singular: 0
  total_forcing_time: 339551 nanoseconds
  total_jacobian_time: 1.60283e+06 nanoseconds
  total_linear_factor_time: 4.88587e+06 nanoseconds
  total_linear_solve_time: 928501 nanoseconds
  JIT result stats: 
  accepted: 162
  function_calls: 334
  jacobian_updates: 162
  number_of_steps: 172
  accepted: 162
  rejected: 0
  decompositions: 172
  solves: 516
  singular: 0
  total_forcing_time: 317043 nanoseconds
  total_jacobian_time: 1.57724e+06 nanoseconds
  total_linear_factor_time: 1.56572e+06 nanoseconds
  total_linear_solve_time: 630749 nanoseconds 