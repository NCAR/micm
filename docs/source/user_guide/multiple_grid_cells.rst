.. _Multiple grid cells:

Multiple Grid Cells
===================

This tutorial will focus on running multiple grid cells. Because the 
:ref:`Rate constants` and :ref:`User defined rate constants` showed both configuring
the mechanism by hand and building with a configuration file, we will only show building
up the mechanism by hand for this tutorial. We will use a simple 3-reaction 3-species mechanism.

.. math::

  A &\longrightarrow B, &k_{1, \mathrm{user\ defined}} \\
  2B &\longrightarrow B + C, &k_{2, \mathrm{user\ defined}} \\
  B + C &\longrightarrow A + C, \qquad &k_{3, \mathrm{user\ defined}} \\

We will use three grid cells. The second grid cells will have concentrations twice as large as the first grid cell.
The third grid cell will have concentrations half as large as the first grid cell. 

If you're looking for a copy and paste, choose
the appropriate tab below and be on your way! Otherwise, stick around for a line by line explanation.

.. tab-set::

    .. tab-item:: Build the Mechanism with the API

      .. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
        :language: cpp
        
Line-by-line explanation
------------------------

This mechanism only needs the user defined rate constant and the rosenbrock solver.

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 1-4

After that, we'll use the ``micm`` namespace.

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 6-7

To create a :cpp:class:`micm::RosenbrockSolver`, we have to define a chemical system (:cpp:class:`micm::System`)
and our reactions, which will be a vector of :cpp:class:`micm::Process` We will use the species to define these as
well as a :cpp:class:`micm::Phase`.

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 11-15


With the species and gas phase, we can define all of our reactions

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 17-33

Now we can define our RosenbrockSolver. This time we'll form the reactions and chemical system in place.
Also, notice the ``false`` in our :cpp:class:`micm::RosenbrockSolverParameters`. This tells the solver
not to reorder the state variables. The reordering is an optimization that can minizie fill-in for some 
of the linear algebra operations. For this example, it's turned off so that the order of the state matches
the order the species are added to the gas phase.

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 35-37


Now we need to get a state and set the concentations of each of the species. In the 
":ref:`Rate constants set concentations`" section of the rate constants tutorial, 
we used a ``std::unordered_map<std::string, std::vector<double>>`` 
to set the concentrations. Here we will set the concentations by providing the :cpp:class:`micm::Species` objects.

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 39-44

Then we set the reaction rates by creating a vector that is 3 elements long, one for each grid cell.

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 46-51

And lastly set the temperature, pressure, and air density for each grid cell.

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 53-62

Now we are ready to run the simulation.

.. literalinclude:: ../../../test/tutorial/test_multiple_grid_cells.cpp
  :language: cpp
  :lines: 64-86


And these are the results.

.. csv-table::
   :header: "time", "grid", "A", "B", "C"
   :widths: 6, 6, 10, 10, 10

        0,     0,   1.00e+00,   0.00e+00,   0.00e+00
        0,     1,   2.00e+00,   0.00e+00,   0.00e+00
        0,     2,   5.00e-01,   0.00e+00,   0.00e+00
      200,     0,   5.35e-01,   4.49e-06,   4.65e-01
      200,     1,   1.21e+00,   6.01e-06,   7.89e-01
      200,     2,   2.30e-01,   3.31e-06,   2.70e-01
      400,     0,   4.50e-01,   3.23e-06,   5.50e-01
      400,     1,   1.05e+00,   4.42e-06,   9.45e-01
      400,     2,   1.85e-01,   2.32e-06,   3.15e-01
      600,     0,   4.00e-01,   2.63e-06,   6.00e-01
      600,     1,   9.59e-01,   3.65e-06,   1.04e+00
      600,     2,   1.60e-01,   1.85e-06,   3.40e-01
      800,     0,   3.64e-01,   2.27e-06,   6.36e-01
      800,     1,   8.90e-01,   3.18e-06,   1.11e+00
      800,     2,   1.42e-01,   1.57e-06,   3.58e-01
      1000,     0,   3.37e-01,   2.01e-06,   6.63e-01
      1000,     1,   8.35e-01,   2.85e-06,   1.16e+00
      1000,     2,   1.29e-01,   1.38e-06,   3.71e-01
      1200,     0,   3.15e-01,   1.82e-06,   6.85e-01
      1200,     1,   7.91e-01,   2.60e-06,   1.21e+00
      1200,     2,   1.19e-01,   1.24e-06,   3.81e-01
      1400,     0,   2.96e-01,   1.67e-06,   7.04e-01
      1400,     1,   7.54e-01,   2.40e-06,   1.25e+00
      1400,     2,   1.11e-01,   1.13e-06,   3.89e-01
      1600,     0,   2.81e-01,   1.55e-06,   7.19e-01
      1600,     1,   7.21e-01,   2.24e-06,   1.28e+00
      1600,     2,   1.04e-01,   1.04e-06,   3.96e-01
      1800,     0,   2.67e-01,   1.45e-06,   7.33e-01
      1800,     1,   6.93e-01,   2.11e-06,   1.31e+00
      1800,     2,   9.77e-02,   9.65e-07,   4.02e-01
      2000,     0,   2.55e-01,   1.37e-06,   7.45e-01
      2000,     1,   6.68e-01,   2.00e-06,   1.33e+00
      2000,     2,   9.25e-02,   9.02e-07,   4.07e-01