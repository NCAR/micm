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

.. tabs::

    .. tab:: Build the Mechanism with the API

      .. literalinclude:: ../../../test/tutorial/test_but_how_fast_is_it.cpp
        :language: cpp
        
Line-by-line explanation
------------------------

Up until now we have neglected to talk about what the solver returns, which is a :cpp:class:`micm::SolverResult`

Now we are ready to run the simulation.

.. literalinclude:: ../../../test/tutorial/test_but_how_fast_is_it.cpp
  :language: cpp
  :lines: 110-133


And these are the results.

.. csv-table::
   :header: "time", "grid", "A", "B", "C"
   :widths: 6, 6, 10, 10, 10

   0, 1, 1.00e+00, 0.00e+00, 0.00e+00
   0, 2, 2.00e+00, 0.00e+00, 0.00e+00
   0, 3, 5.00e-01, 0.00e+00, 0.00e+00
   200, 1, 5.35e-01, 4.49e-06, 4.65e-01
   200, 2, 1.21e+00, 6.01e-06, 7.89e-01
   200, 3, 2.30e-01, 3.31e-06, 2.70e-01
   400, 1, 4.50e-01, 3.23e-06, 5.50e-01
   400, 2, 1.05e+00, 4.42e-06, 9.45e-01
   400, 3, 1.85e-01, 2.32e-06, 3.15e-01
   600, 1, 4.00e-01, 2.63e-06, 6.00e-01
   600, 2, 9.59e-01, 3.65e-06, 1.04e+00
   600, 3, 1.60e-01, 1.85e-06, 3.40e-01
   800, 1, 3.64e-01, 2.27e-06, 6.36e-01
   800, 2, 8.90e-01, 3.18e-06, 1.11e+00
   800, 3, 1.42e-01, 1.57e-06, 3.58e-01
   1000, 1, 3.37e-01, 2.01e-06, 6.63e-01
   1000, 2, 8.35e-01, 2.85e-06, 1.16e+00
   1000, 3, 1.29e-01, 1.38e-06, 3.71e-01
   1200, 1, 3.15e-01, 1.82e-06, 6.85e-01
   1200, 2, 7.91e-01, 2.60e-06, 1.21e+00
   1200, 3, 1.19e-01, 1.24e-06, 3.81e-01
   1400, 1, 2.96e-01, 1.67e-06, 7.04e-01
   1400, 2, 7.54e-01, 2.40e-06, 1.25e+00
   1400, 3, 1.11e-01, 1.13e-06, 3.89e-01
   1600, 1, 2.81e-01, 1.55e-06, 7.19e-01
   1600, 2, 7.21e-01, 2.24e-06, 1.28e+00
   1600, 3, 1.04e-01, 1.04e-06, 3.96e-01
   1800, 1, 2.67e-01, 1.45e-06, 7.33e-01
   1800, 2, 6.93e-01, 2.11e-06, 1.31e+00
   1800, 3, 9.77e-02, 9.65e-07, 4.02e-01
   2000, 1, 2.55e-01, 1.37e-06, 7.45e-01
   2000, 2, 6.68e-01, 2.00e-06, 1.33e+00
   2000, 3, 9.25e-02, 9.02e-07, 4.07e-01