.. _Multiple grid cells:

Multiple Grid Cells
===================

This tutorial will focus on running multiple grid cells. Because the 
:ref:`Rate constants` and :ref:`User defined rate constants` showed both configuring
the mechanism by hand and building with a configuration file, we will only show building
up the mechanism by hand for this tutorial. That mechanism will be a simple chapman mechanism.


.. math::

  O_{1}^{d} + N_2 &\longrightarrow O + N_2, &k_{1, \mathrm{arrhenius}} \\
  O_{1}^{d} + O_2 &\longrightarrow O + O_2, &k_{2, \mathrm{arrhenius}} \\
  O + O_3 &\longrightarrow 2O_2, &k_{3, \mathrm{arrhenius}} \\
  O + O_2 + M &\longrightarrow O_3 + M, &k_{4, \mathrm{arrhenius}} \\
  O_2 &\longrightarrow 2O, &k_{5, \mathrm{photolysis}} \\
  O_3 &\longrightarrow O_1^d +O_2,\qquad &k_{6, \mathrm{photolysis}} \\
  O_3 &\longrightarrow O +O_2, &k_{7, \mathrm{photolysis}} \\

We will use three grid cells. The second grid cells will have concentrations twice as large as the first grid cell.
The third grid cell will have concentrations half as large as the first grid cell.