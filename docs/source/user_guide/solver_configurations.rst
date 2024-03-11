.. _Solver configurations:

Solver configurations
=====================

This tutorial will focus on configuring the solver with different stages.
We will use a simple 3-reaction 3-species mechanism. The setup here is the same in
:ref:`Multiple grid cells`, except only one grid cell is used.

.. math::

  A &\longrightarrow B, &k_{1, \mathrm{user\ defined}} \\
  2B &\longrightarrow B + C, &k_{2, \mathrm{user\ defined}} \\
  B + C &\longrightarrow A + C, \qquad &k_{3, \mathrm{user\ defined}} \\

If you're looking for a copy and paste, choose
the appropriate tab below and be on your way! Otherwise, stick around for a line by line explanation.

.. tab-set::

    .. tab-item:: Build the Mechanism with the API

      .. literalinclude:: ../../../test/tutorial/test_solver_configuration.cpp
        :language: cpp
        
Line-by-line explanation
------------------------

There are a total of five different configurations for the Rosenbrock solver. Each corresponds to a different number of stages.
What each configuration means is beyond the scope of this tutorial. However, references are included for further reading.

- Two stage
    - An L-stable method, 2 stages, order 2. :cite:`Hairer1996`
- Three stage
    - An L-stable method, 3 stages, order 3. :cite:`SanduStiff2-1997`
- Four stage
    - An L-stable method, 4 stages, order 4. :cite:`Hairer1996`
- Four stage, differential algebraic
    - A stiffly stable method, 4 stages, order 3, :cite:`Hairer1996`
- Six stage, differential algebraic
    - Stiffly stable rosenbrock method of order 4 with 6 stages, :cite:`Hairer1996`


Configuring the rosenbrock solver is as easy as providing the solver with a set of :cpp:class:`micm::RosenbrockSolverParameters`

.. literalinclude:: ../../../test/tutorial/test_solver_configuration.cpp
  :language: cpp
  :lines: 121-133

After that, the usage is the same as a regular solver. A templated method was used here to run the same solving code
for each of the different solver configurations.

.. literalinclude:: ../../../test/tutorial/test_solver_configuration.cpp
  :language: cpp
  :lines: 10-90

Running this program should give an output similar to this:

.. code-block:: console

  Two stages: 
  time,         A,         B,         C
      0,  1.00e+00,  0.00e+00,  0.00e+00
    200,  5.37e-01,  4.50e-06,  4.63e-01
    400,  4.51e-01,  3.23e-06,  5.49e-01
    600,  4.01e-01,  2.64e-06,  5.99e-01
    800,  3.65e-01,  2.28e-06,  6.35e-01
  1000,  3.38e-01,  2.02e-06,  6.62e-01
  1200,  3.16e-01,  1.83e-06,  6.84e-01
  1400,  2.97e-01,  1.68e-06,  7.03e-01
  1600,  2.82e-01,  1.56e-06,  7.18e-01
  1800,  2.68e-01,  1.46e-06,  7.32e-01
  2000,  2.56e-01,  1.37e-06,  7.44e-01
  Total solve time: 135748 nanoseconds
  accepted: 178
  function_calls: 356
  jacobian_updates: 178
  number_of_steps: 178
  accepted: 178
  rejected: 0
  decompositions: 178
  solves: 356
  singular: 0
  total_forcing_time: 9962 nanoseconds
  total_jacobian_time: 6454 nanoseconds
  total_linear_factor_time: 26835 nanoseconds
  total_linear_solve_time: 12044 nanoseconds

  Three stages: 
  time,         A,         B,         C
      0,  1.00e+00,  0.00e+00,  0.00e+00
    200,  5.35e-01,  4.50e-06,  4.65e-01
    400,  4.50e-01,  3.23e-06,  5.50e-01
    600,  4.00e-01,  2.63e-06,  6.00e-01
    800,  3.64e-01,  2.27e-06,  6.36e-01
  1000,  3.37e-01,  2.01e-06,  6.63e-01
  1200,  3.15e-01,  1.82e-06,  6.85e-01
  1400,  2.96e-01,  1.67e-06,  7.04e-01
  1600,  2.81e-01,  1.55e-06,  7.19e-01
  1800,  2.67e-01,  1.45e-06,  7.33e-01
  2000,  2.55e-01,  1.37e-06,  7.45e-01
  Total solve time: 110002 nanoseconds
  accepted: 127
  function_calls: 254
  jacobian_updates: 127
  number_of_steps: 127
  accepted: 127
  rejected: 0
  decompositions: 127
  solves: 381
  singular: 0
  total_forcing_time: 7421 nanoseconds
  total_jacobian_time: 4295 nanoseconds
  total_linear_factor_time: 19248 nanoseconds
  total_linear_solve_time: 12093 nanoseconds

  Four stages: 
  time,         A,         B,         C
      0,  1.00e+00,  0.00e+00,  0.00e+00
    200,  5.36e-01,  4.49e-06,  4.64e-01
    400,  4.50e-01,  3.20e-06,  5.50e-01
    600,  4.00e-01,  2.62e-06,  6.00e-01
    800,  3.64e-01,  2.26e-06,  6.36e-01
  1000,  3.37e-01,  2.01e-06,  6.63e-01
  1200,  3.15e-01,  1.82e-06,  6.85e-01
  1400,  2.97e-01,  1.67e-06,  7.03e-01
  1600,  2.81e-01,  1.55e-06,  7.19e-01
  1800,  2.67e-01,  1.45e-06,  7.33e-01
  2000,  2.56e-01,  1.37e-06,  7.44e-01
  Total solve time: 127291 nanoseconds
  accepted: 127
  function_calls: 381
  jacobian_updates: 127
  number_of_steps: 127
  accepted: 127
  rejected: 0
  decompositions: 127
  solves: 508
  singular: 0
  total_forcing_time: 9749 nanoseconds
  total_jacobian_time: 4537 nanoseconds
  total_linear_factor_time: 20040 nanoseconds
  total_linear_solve_time: 17923 nanoseconds

  Four stages differential algebraic: 
  time,         A,         B,         C
      0,  1.00e+00,  0.00e+00,  0.00e+00
    200,  5.36e-01,  4.49e-06,  4.64e-01
    400,  4.51e-01,  3.23e-06,  5.49e-01
    600,  4.00e-01,  2.63e-06,  6.00e-01
    800,  3.64e-01,  2.27e-06,  6.36e-01
  1000,  3.37e-01,  2.01e-06,  6.63e-01
  1200,  3.15e-01,  1.82e-06,  6.85e-01
  1400,  2.97e-01,  1.68e-06,  7.03e-01
  1600,  2.81e-01,  1.55e-06,  7.19e-01
  1800,  2.68e-01,  1.45e-06,  7.32e-01
  2000,  2.56e-01,  1.37e-06,  7.44e-01
  Total solve time: 126334 nanoseconds
  accepted: 128
  function_calls: 384
  jacobian_updates: 128
  number_of_steps: 128
  accepted: 128
  rejected: 0
  decompositions: 128
  solves: 512
  singular: 0
  total_forcing_time: 10584 nanoseconds
  total_jacobian_time: 4376 nanoseconds
  total_linear_factor_time: 19792 nanoseconds
  total_linear_solve_time: 17254 nanoseconds


  Six stages differential algebraic: 
    time,         A,         B,         C
        0,  1.00e+00,  0.00e+00,  0.00e+00
      200,  5.36e-01,  4.49e-06,  4.64e-01
      400,  4.51e-01,  3.22e-06,  5.49e-01
      600,  4.00e-01,  2.63e-06,  6.00e-01
      800,  3.64e-01,  2.27e-06,  6.36e-01
    1000,  3.37e-01,  2.01e-06,  6.63e-01
    1200,  3.15e-01,  1.82e-06,  6.85e-01
    1400,  2.97e-01,  1.67e-06,  7.03e-01
    1600,  2.81e-01,  1.55e-06,  7.19e-01
    1800,  2.67e-01,  1.45e-06,  7.33e-01
    2000,  2.56e-01,  1.37e-06,  7.44e-01
  Total solve time: 174959 nanoseconds
  accepted: 127
  function_calls: 762
  jacobian_updates: 127
  number_of_steps: 127
  accepted: 127
  rejected: 0
  decompositions: 127
  solves: 762
  singular: 0
  total_forcing_time: 20743 nanoseconds
  total_jacobian_time: 4497 nanoseconds
  total_linear_factor_time: 20085 nanoseconds
  total_linear_solve_time: 28233 nanoseconds