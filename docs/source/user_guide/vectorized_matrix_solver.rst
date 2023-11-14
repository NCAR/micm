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

Up until now, we've mostly created our :cpp:class:`micm::RosenbrockSolver` solvers like this:

.. code-block:: cpp
   
     RosenbrockSolver<> solver{ System(SystemParameters{ .gas_phase_ = gas_phase }),
                             std::vector<Process>{ r1, r2, r3 },
                             RosenbrockSolverParameters::three_stage_rosenbrock_parameters(3, false) };

The empty angle brackets ``<>`` tell the compiler to use the default template paramters. The first two are most important and
denote the ``MatrixPolicy`` and ``SparseMatrixPolicy`` to be used when we solve. By default, these are :cpp:class:`micm::Matrix`
and :cpp:class:`micm::SparseMatrix`. The most visible data affected by these paratemeters is the state's 
:cpp:member:`micm::State::variables_` parameter, representing the concnetration of each of the species. 