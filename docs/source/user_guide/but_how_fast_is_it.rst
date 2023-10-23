.. _But how fast is it:

But how fast is it?
===================

This tutorial will focus on timing the solver to show how you can measure performance.
We will use a simple 3-reaction 3-species mechanism. The setup here is the same in
:ref:`Multiple grid cells`. To understand the full setup, read that tutorial. Otherwise,
we assume that configuring a rosenbrock solver is understood and instead we will focus on timing 
the solver.

.. math::

  A &\longrightarrow B, &k_{1, \mathrm{user\ defined}} \\
  2B &\longrightarrow B + C, &k_{2, \mathrm{user\ defined}} \\
  B + C &\longrightarrow A + C, \qquad &k_{3, \mathrm{user\ defined}} \\

If you're looking for a copy and paste, choose
the appropriate tab below and be on your way! Otherwise, stick around for a line by line explanation.

.. tab-set::

    .. tab-item:: Build the Mechanism with the API

      .. literalinclude:: ../../../test/tutorial/test_but_how_fast_is_it.cpp
        :language: cpp
        
Line-by-line explanation
------------------------

Up until now we have neglected to talk about what the solver returns, which is a :cpp:class:`micm::RosenbrockSolver::SolverResult`.

There are four values returned.

#. :cpp:member:`micm::RosenbrockSolver::SolverResult::final_time_`

    * This is the final simulation time achieved by the solver. The :cpp:func:`micm::RosenbrockSolver::Solve` function attempts to integrate the passed in state forward a set number of seconds. Often, the solver is able to complete the integration.  However, extremely stiff systems may only solve for a fraction of the time. It is imperative that the ``final_time_`` value is checked. If it is not equal to the amount of time you intended to solve for, call solve again as we do in the tutorials with the difference between what was solved and how long you intended to solve.

      .. note::
        This does **not** represent the amount of time taken by the solve routine. You must measure that yourself(shown below). The ``final_time_`` is simulation time.

#. :cpp:member:`micm::RosenbrockSolver::SolverResult::result_`

    * This contains the integrated state value; the concentrations reached at the end of Solve function after the amount of time specified by ``final_time_``.

#. :cpp:member:`micm::RosenbrockSolver::SolverResult::state_`

    * There are many possible reasons for the solver to return. This value is one of the possible enum values define on the :cpp:enum:`micm::SolverState`. Hopefully, you receive a :cpp:enumerator:`micm::SolverState::Converged` state. But, it is good to always check this to ensure the solver really did converge. You can print this value using the :cpp:func:`micm::StateToString` function.

#. :cpp:member:`micm::RosenbrockSolver::SolverResult::stats_`

    * This is an instance of a :cpp:class:`micm::RosenbrockSolver::SolverStats` struct which contains information about the number of function calls and optionally the total cumulative amount of time spent calling each function. For the time to be collected, you must call the ``Solve`` function with a ``true`` templated parameter. Please see the example below.

First, let's run the simulation but without collecting the solve time. We'll inspect the solver state and look at what's collected
in the stats object. 

.. literalinclude:: ../../../test/tutorial/test_but_how_fast_is_it.cpp
  :language: cpp
  :lines: 69-80

.. code-block:: bash

  Solver state: Converged
  accepted: 20
  function_calls: 40
  jacobian_updates: 20
  number_of_steps: 20
  accepted: 20
  rejected: 0
  decompositions: 20
  solves: 60
  singular: 0

To get the total accumulated time of each function call, you need to specify the templated boolean argument to turn the timing on.
We can also record the total runtime of the ``Solve`` function. Through the magic of templates, the timing information is only
collected when you use the ``true`` version of the templated function. These values are not even computed, meaning no CPU cycles are 
wasted, for the ``false`` version.

.. literalinclude:: ../../../test/tutorial/test_but_how_fast_is_it.cpp
  :language: cpp
  :lines: 82-90

.. code-block:: bash

  Total solve time: 24416 nanoseconds
  total_forcing_time: 3167 nanoseconds
  total_jacobian_time: 1710 nanoseconds
  total_linear_factor_time: 4584 nanoseconds
  total_linear_solve_time: 3290 nanoseconds

.. note::
  Your systems clock may not report the same granularity depending on how accurate your system clock is. 
