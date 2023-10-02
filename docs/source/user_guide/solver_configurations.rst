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

.. tabs::

    .. tab:: Build the Mechanism with the API

      .. literalinclude:: ../../../test/tutorial/test_solver_configuration.cpp
        :language: cpp
        
Line-by-line explanation
------------------------

There are a total of six different configurations for the Rosenbrock solver. Each corresponds to a different number of stages.
What each configuration means is beyond the scope of this tutorial. However, more information can be found in .

:cite:`Sandu1996`

.. literalinclude:: ../../../test/tutorial/test_solver_configuration.cpp
  :language: cpp
  :lines: 88-96