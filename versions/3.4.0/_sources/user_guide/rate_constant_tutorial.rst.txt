.. _Rate constants:

Rate Constants (except user defined ones)
=========================================

MICM supports a subset of the rate constants defined as part of the 
`OpenAtmos Mechanism Configuration <https://open-atmos.github.io/MechanismConfiguration/reactions/index.html>`_
We will be adding more in the future. 
The links to the ``micm`` classes below detail the class and methods. Please check the OpenAtmos standard for
specific on the algorithm implemented for each rate constant type.
At present, supported rate constants are:

- :cpp:class:`micm::ArrheniusRateConstant`, `OpenAtmos Arrhenius description <https://open-atmos.github.io/MechanismConfiguration/reactions/arrhenius.html>`_
- :cpp:class:`micm::BranchedRateConstant`, `OpenAtmos Branched description <https://open-atmos.github.io/MechanismConfiguration/reactions/branched.html>`_
- :cpp:class:`micm::SurfaceRateConstant`, `OpenAtmos Surface description <https://open-atmos.github.io/MechanismConfiguration/reactions/surface.html>`_
- :cpp:class:`micm::TernaryChemicalActivationRateConstant`, `OpenAtmos Ternay chemical activiation description <https://open-atmos.github.io/MechanismConfiguration/reactions/ternary_chemical_activation.html>`_
- :cpp:class:`micm::TroeRateConstant`, `OpenAtmos Troe description <https://open-atmos.github.io/MechanismConfiguration/reactions/troe.html>`_
- :cpp:class:`micm::TunnelingRateConstant`, `OpenAtmos Tunneling description <https://open-atmos.github.io/MechanismConfiguration/reactions/tunneling.html>`_
- :cpp:class:`micm::UserDefinedRateConstant`, this is a custom type, but this is how photolysis rates are setup `OpenAtmos Photolysis description <https://open-atmos.github.io/MechanismConfiguration/reactions/photolysis.html>`_

This tutorial covers all but the last one. See the :ref:`User defined rate constants` tutorial for examples and use
cases on that.

We'll setup and solve a fake chemical system with 7 species and 6 reactions, 

.. math::

  A &\longrightarrow B, &k_{1, \mathrm{arrhenius}} \\
  B &\longrightarrow C (\mathrm{alkoxy\ products}) + D (\mathrm{nitrate\ products}), &k_{2, \mathrm{branched}} \\
  C &\longrightarrow E, &k_{3, \mathrm{surface}} \\
  D &\longrightarrow 2F, &k_{4, \mathrm{ternary\ chemical\ activation}} \\
  2E &\longrightarrow G, &k_{5, \mathrm{troe}} \\
  F &\longrightarrow G, &k_{6, \mathrm{tunneling}} \\


MICM can be configured in two ways. We can either build the mechanism up by hand with the ``micm`` API,
or parse a valid mechanism Configuration
in the OpenAtmos format. In this tutorial, we will do both. 

If you're looking for a copy and paste, choose
the appropriate tab below and be on your way! Otherwise, stick around for a line by line explanation.


.. tab-set::

    .. tab-item:: Build the Mechanism with the API

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp

    .. tab-item:: OpenAtmos Configuration reading

        .. raw:: html

            <div class="download-div">
              <a href="../_static/tutorials/rate_constants_no_user_defined.zip" download>
                <button class="download-button">Download zip configuration</button>
              </a>
            </div>

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_with_config.cpp
          :language: cpp

Line-by-line explanation
------------------------

To get started, we'll need to include each rate constant type and the 
rosenbrock solver at the top of the file.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 1-13

After that, we'll use the ``micm`` namespace so that we don't have to repeat it everywhere we need it.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 14-15

To create a :cpp:class:`micm::RosenbrockSolver`, we have to define a chemical system (:cpp:class:`micm::System`)
and our reactions, which will be a vector of :cpp:class:`micm::Process` We will use the species to define these.

.. tab-set::

    .. tab-item:: Build the Mechanism with the API

        To do this by hand, we have to define all of the chemical species in the system. This allows us to set
        any properties of the species that may be necessary for rate constanta calculations, like molecular weights 
        and diffusion coefficients for the surface reaction.  We will also put these species into the gas phase.

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 19-30

        Now that we have a gas phase and our species, we can start building the reactions. Two things to note are that
        stoichiemtric coefficients for reactants are represented by repeating that product as many times as you need.
        To specify the yield of a product, we've created a typedef :cpp:type:`micm::Yield` 
        and a function :cpp:func:`micm::yields` that produces these.

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 32-96
        
        And finally we define our chemical system and reactions

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 98-99

    .. tab-item:: OpenAtmos Configuration reading

        After defining a valid OpenAtmos configuration with reactions that ``micm`` supports, configuring the chemical
        system and the processes is as simple as using the :cpp:class:`micm::SolverConfig` class

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_with_config.cpp
          :language: cpp
          :lines: 20-32

Now that we have a chemical system and a list of reactions, we can create the RosenbrockSolver.
There are several ways to configure the solver. Here we are using a three stage solver. More options
can be found in the :cpp:class:`micm::RosenbrockSolverParameters` and in the :ref:`Solver Configurations` tutorial.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 101

The rosenbrock solver will provide us a state, which we can use to set the concentrations,
custom rate parameters, and temperature and pressure. Note that setting the custom rate paramters is different depending
on if you define the configuration by hand or read it in. The parser has defaults for the names of the custom parameters
and when defined by hand we choose these.

.. _Rate constants set concentations:

Initializing the state
^^^^^^^^^^^^^^^^^^^^^^

.. tab-set::

    .. tab-item:: Build the Mechanism with the API

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 102-116

    .. tab-item:: OpenAtmos Configuration reading

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_with_config.cpp
          :language: cpp
          :lines: 36-54


Finally, we are ready to pick a timestep and solve the system.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 118-142


This is the output:


.. csv-table:: Table Title
   :header: "time", "A", "B", "C", "D", "E", "F", "G"
   :widths: 10, 15, 15, 15, 15, 15, 15, 15

      0,   1.00e+00,   0.00e+00,   0.00e+00,   0.00e+00,   0.00e+00,   0.00e+00,   0.00e+00
    500,   3.18e-09,   3.66e-09,   9.83e-01,   3.88e-14,   1.41e-03,   2.02e-13,   7.92e-03
    1000,   1.14e-14,   1.31e-14,   9.66e-01,   1.39e-19,   1.40e-03,   7.24e-19,   1.64e-02
    1500,   4.09e-20,   4.71e-20,   9.49e-01,   4.98e-25,   1.39e-03,   2.59e-24,   2.48e-02
    2000,   1.47e-25,   1.69e-25,   9.33e-01,   1.79e-30,   1.38e-03,   9.30e-30,   3.30e-02
    2500,   5.26e-31,   6.05e-31,   9.17e-01,   6.40e-36,   1.37e-03,   3.33e-35,   4.11e-02
    3000,   1.89e-36,   2.17e-36,   9.01e-01,   2.30e-41,   1.36e-03,   1.20e-40,   4.90e-02
    3500,   6.77e-42,   7.78e-42,   8.85e-01,   8.23e-47,   1.34e-03,   4.29e-46,   5.68e-02
    4000,   2.43e-47,   2.79e-47,   8.70e-01,   2.95e-52,   1.33e-03,   1.54e-51,   6.44e-02
    4500,   8.70e-53,   1.00e-52,   8.55e-01,   1.06e-57,   1.32e-03,   5.51e-57,   7.20e-02
    5000,   3.12e-58,   3.59e-58,   8.40e-01,   3.80e-63,   1.31e-03,   1.98e-62,   7.94e-02
