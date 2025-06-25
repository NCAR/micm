.. _OpenMP:

OpenMP
======

This tutorial will focus on running micm with OpenMP support.
We will use a simple 3-reaction 3-species mechanism. The setup here is the same in
:ref:`Multiple grid cells`, except only one grid cell is used.

.. math::
  A &\longrightarrow B, &k_{1, \mathrm{user\ defined}} \\
  2B &\longrightarrow B + C, &k_{2, \mathrm{user\ defined}} \\
  B + C &\longrightarrow A + C, \qquad &k_{3, \mathrm{user\ defined}} \\

If you're looking for a copy and paste, copy below and be on your way! Otherwise, stick around for a line by line explanation.

.. tab-set::

    .. tab-item:: OpenAtmos Configuration reading

      .. raw:: html

          <div class="download-div">
            <a href="../_static/tutorials/robertson.zip" download>
              <button class="download-button">Download zip configuration</button>
            </a>
          </div>

      .. literalinclude:: ../../../test/tutorial/test_openmp.cpp
        :language: cpp

Line-by-line explanation
------------------------

At present, the Rosenbrock class is not thread safe. Thread safety is being added and will be available soon in new release.

Until then, you can run one instance of the Rosenbrock solver on its own thread. Configuration data, at least, can be shared
across the threads. This tutorial reads a configuraiton file and then configures three solvers, each on their own thread,
with that same configuration information.

First we need to bring in the necessary imports and use the micm namespace.

.. literalinclude:: ../../../test/tutorial/test_openmp.cpp
  :language: cpp
  :lines: 1-6

Then we'll define a funtion that we can call with the configuration information from each thread. This function will build
and run the rosenbrock solver.

.. literalinclude:: ../../../test/tutorial/test_openmp.cpp
  :language: cpp
  :lines: 26-63

The main function simply reads the configuration file, sets the number of threads, and then sets up an OpenMP blcok
with three threads. The function defined above is called on each thread.

.. literalinclude:: ../../../test/tutorial/test_openmp.cpp
  :language: cpp
  :lines: 65-107

Running this program should give an output similar to this:

.. code-block:: console

  Thread 1
           A,         B,         C
    2.55e-01,  1.37e-06,  7.45e-01
  Thread 2
           A,         B,         C
    2.55e-01,  1.37e-06,  7.45e-01
  Thread 3
           A,         B,         C
    2.55e-01,  1.37e-06,  7.45e-01