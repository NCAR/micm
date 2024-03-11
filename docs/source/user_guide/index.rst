##########
User Guide
##########

MICM is quite expansive in its capabilities. We've written examples showing many different use cases.
Hopefully our examples are exhaustive enough to provide you with the tools you need to solve some atmospheric chemistry.

If you happen to find our examples are lacking for your needs, please, 
`fill out an issue <https://github.com/NCAR/micm/issues/new>`_ and request the kind of example you'd like. 


All of these tutorials are included in our automated tests. Each of them can be found in the code base in the
``test/tutorial`` directory. When building MICM with tests (the default), you can each test individually to see the output.


.. code-block:: console
  
    $ git clone https://github.com/NCAR/micm.git
    $ cd micm
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ ./test_multiple_grid_cells
    $ ./test_rate_constants_no_user_defined_example_by_hand
    $ ./test_rate_constants_user_defined_example_by_hand


If you would like to include the json examples, you must configure micm to build with json support.

.. code-block:: console
  
    $ cmake -DMICM_ENABLE_JSON=ON ..
    $ make
    $ ./test_rate_constants_no_user_defined_example_with_config
    $ ./test_rate_constants_user_defined_example_with_config


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation_and_usage
   rate_constant_tutorial
   user_defined_rate_constant_tutorial
   multiple_grid_cells
   solver_configurations
   but_how_fast_is_it
   openmp
   jit
   vectorized_matrix_solver
