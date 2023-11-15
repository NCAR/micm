.. _Vectorized matrix solver:

Vectorized matrix solver
########################

Many of the tutorials have defined a chemical system without talking about their location. Each of these is essentially
defining a `box model <https://en.wikipedia.org/wiki/Climate_model#Box_models>`_. We often want to model chemistry
in more than one location; multiple boxes. In weather and climate models, the domain is broken up into grids 2d or 3d grids.
Each location may be referred to as a grid cell or column. With ``micm``, we can solve multiple grid cells simultaneously in
several ways.


The tutorials listed below solve chemistry in more than one grid cell.

* :ref:`Multiple grid cells`, this tutorial does so with a different data structure
* :ref:`OpenMP`, this tutorial uses a single thread per grid cell


In :ref:`Multiple grid cells` we assumed we were solving for the concentrations of an :math:`N\times D` matrix, :math:`N` species and 
:math:`D` grid cells. However, we used the default matrix ordering in that tutorial. This is fine for simulations that aren't
too large, but computationally may be inferior. This is because the same species appears in each grid cell, :math:`N` elements
away from each other in memory. When solving the chemical equations, the same mathematical operations have to be applied to each
species. We can organize the memory of the matrix so that the same species across all :math:`D` grid cells are directly next
to each other in memory. In this way, we help to encourage the compiler to vectorize the code and take advantage of 
`SIMD instructions <https://en.wikipedia.org/wiki/Single_instruction,_multiple_data>`_.

Organizing the data this way may be familiar to GPU programmers as
`coalesced memory access <https://www.sciencedirect.com/topics/computer-science/memory-access-coalescing>`_.

To understand this, we'll look at how the memeory is organized when solving the robertson equations with 3 grid cells using
the various matrix types and ordering in ``micm``.

.. math::

  A &\longrightarrow B, &k_{1, \mathrm{user\ defined}} \\
  2B &\longrightarrow B + C, &k_{2, \mathrm{user\ defined}} \\
  B + C &\longrightarrow A + C, \qquad &k_{3, \mathrm{user\ defined}} \\


This entire example is contained below. Copy and paste to follow along.

.. tab-set::

    .. tab-item:: Build the Mechanism with the API

        .. literalinclude:: ../../../test/tutorial/test_vectorized_matrix_solver.cpp
          :language: cpp


Matrix types in micm
--------------------

``micm`` has its own matrix types, each providing a separate data layout. The data layouts are meant to provide efficient access
such that cache misses are minimized. This is done by placing all of the data into a contiguous vector, and optionally imposing
some sort of ordering on the data, like grouping species across gridcells together in memory. 


The dimensions of all ``micm``
matrices are :math:`\mathrm{number\ of\ grid\ cells} \times \mathrm{number\ of\ species}`

* :cpp:class:`micm::Matrix`, a dense matrix, what you probably think of when you think of a matrix, with the caveat that all rows and columns live in one long array

* :cpp:class:`micm::VectorMatrix`, the matrix type of interest here, where the data is ordered such that species across grid cells lie beside each other in memory

* :cpp:class:`micm::SparseMatrix`, a `sparse matrix <https://en.wikipedia.org/wiki/Sparse_matrix>`_ data structure

   * The sparse matrix also allows you to specify wether the data is ordered as a regular dense matrix, with zeros removed, or with a vectorized ordering, the topic of this tutorial.

      * :cpp:class:`micm::SparseMatrixStandardOrdering`

      * :cpp:class:`micm::SparseMatrixVectorOrdering`

Regular, dense matrix ordering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: images/dense.png

Up until now, we've mostly created our :cpp:class:`micm::RosenbrockSolver` solvers like this:

.. literalinclude:: ../../../test/tutorial/test_vectorized_matrix_solver.cpp
  :language: cpp
  :lines: 86

The empty angle brackets ``<>`` tell the compiler to use the default template paramters. The first two are most important and
denote the ``MatrixPolicy`` and ``SparseMatrixPolicy`` to be used when we solve. By default, these are :cpp:class:`micm::Matrix`
and :cpp:class:`micm::SparseMatrix`. The most visible data affected by these paratemeters is the state's 
:cpp:member:`micm::State::variables_` parameter, representing the concnetration of each of the species. 

After we created that first ``solver``, if we get a state and set the concentration of the species and then print them,
they are printed in order of grid cell first and then species. Here, when setting the species, read the double as ``grid cell.species``.

.. literalinclude:: ../../../test/tutorial/test_vectorized_matrix_solver.cpp
  :language: cpp
  :lines: 88-97

.. code-block:: bash

  1.1
  1.2
  1.3
  2.1
  2.2
  2.3
  3.1
  3.2
  3.3

Vectorized matrix ordering
^^^^^^^^^^^^^^^^^^^^^^^^^^

But, if we create a vectorized solver, set the same concentrations and print out the state variables, we see that it is organized
first by species and then by grid cell. Note that we needed to set some partial template speciailizations at the top of the file
to create these.

At the top of the file:

.. literalinclude:: ../../../test/tutorial/test_vectorized_matrix_solver.cpp
  :language: cpp
  :lines: 14-18

and then we pass these templates as the template arguments to the vectorized Rosenbrock solver 

.. literalinclude:: ../../../test/tutorial/test_vectorized_matrix_solver.cpp
  :language: cpp
  :lines: 101-112

.. code-block:: bash

  1.1
  2.1
  3.1
  1.2
  2.2
  3.2
  1.3
  2.3
  3.3

And that's all you have to do to orgnaize the data by species first. By specifying the template parameter on the 
solver, each operation will use the same ordering for all of the data structures needed to solve the chemical system.


Sparse matrix, standard ordering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: images/sparse.png

Sparse matrix, vector ordering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use the vectorized just as you would the regular solver, and in fact the same output is produced.

.. literalinclude:: ../../../test/tutorial/test_vectorized_matrix_solver.cpp
  :language: cpp
  :lines: 116-134

.. code-block:: bash

  Cell 0
   Species             Regular          Vectorized
         A             1.02174             1.02174
         B         1.58226e-06         1.58226e-06
         C             2.57826             2.57826

  Cell 1
    Species             Regular          Vectorized
          A             2.01153             2.01153
          B         1.75155e-06         1.75155e-06
          C             4.58846             4.58846

  Cell 2
    Species             Regular          Vectorized
          A              3.0096              3.0096
          B         1.82514e-06         1.82514e-06
          C              6.5904              6.5904

