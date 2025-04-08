.. _Rate constants:

Rate Constants (except user defined ones)
=========================================

MICM supports a subset of the rate constants defined as part of the 
`OpenAtmos Mechanism Configuration <https://open-atmos.github.io/MechanismConfiguration/reactions/index.html>`_
We will be adding more in the future. 
The links to the ``micm`` classes below detail the class and methods. Please check the OpenAtmos standard for
specifics on the algorithm implemented for each rate constant type.
At present, supported rate constants are:

- :cpp:class:`micm::ArrheniusRateConstant`, `OpenAtmos Arrhenius description <https://open-atmos.github.io/MechanismConfiguration/reactions/arrhenius.html>`_
- :cpp:class:`micm::BranchedRateConstant`, `OpenAtmos Branched description <https://open-atmos.github.io/MechanismConfiguration/reactions/branched.html>`_
- :cpp:class:`micm::SurfaceRateConstant`, `OpenAtmos Surface description <https://open-atmos.github.io/MechanismConfiguration/reactions/surface.html>`_
- :cpp:class:`micm::TernaryChemicalActivationRateConstant`, `OpenAtmos Ternay chemical activiation description <https://open-atmos.github.io/MechanismConfiguration/reactions/ternary_chemical_activation.html>`_
- :cpp:class:`micm::TroeRateConstant`, `OpenAtmos Troe description <https://open-atmos.github.io/MechanismConfiguration/reactions/troe.html>`_
- :cpp:class:`micm::TunnelingRateConstant`, `OpenAtmos Tunneling description <https://open-atmos.github.io/MechanismConfiguration/reactions/tunneling.html>`_
- :cpp:class:`micm::UserDefinedRateConstant`, this is a custom type, but this is how photolysis, first order loss, and emission rates are setup `OpenAtmos Photolysis description <https://open-atmos.github.io/MechanismConfiguration/reactions/photolysis.html>`_

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

Line-by-line explanation
------------------------

To get started, we'll need to include each rate constant type and the 
rosenbrock solver at the top of the file.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 1-14

After that, we'll use the ``micm`` namespace so that we don't have to repeat it everywhere we need it.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 16-17

To create a :cpp:class:`micm::RosenbrockSolver`, we have to define a chemical system (:cpp:class:`micm::System`)
and our reactions, which will be a vector of :cpp:class:`micm::Process` We will use the species to define these.

.. tab-set::

    .. tab-item:: Build the Mechanism with the API

        To do this by hand, we have to define all of the chemical species in the system. This allows us to set
        any properties of the species that may be necessary for rate constanta calculations, like molecular weights 
        and diffusion coefficients for the surface reaction.  We will also put these species into the gas phase.

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 25-36

        Now that we have a gas phase and our species, we can start building the reactions. Two things to note are that
        stoichiemtric coefficients for reactants are represented by repeating that product as many times as you need.
        To specify the yield of a product, we've created a typedef :cpp:type:`micm::Yield` 
        and a function :cpp:func:`micm::Yields` that produces these. Note that we add a conversion for
        some rate constant parameters to be consistent with the configuration file that expects rate constants
        to be in cm^3/molecule/s. (All units will be mks in the next version of the configuration file format.)

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 38-102
        
        And finally we define our chemical system and reactions

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 104-105

Now that we have a chemical system and a list of reactions, we can create the RosenbrockSolver.
There are several ways to configure the solver. Here we are using a three stage solver. More options
can be found in the :cpp:class:`micm::RosenbrockSolverParameters` and in the :ref:`Solver Configurations` tutorial.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 107-110

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
          :lines: 111-126


Finally, we are ready to pick a timestep and solve the system.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 129-151


This is the output:


.. csv-table:: The Change of Concentration with Time
   :header: "time", "A", "B", "C", "D", "E", "F", "G"
   :widths: 10, 15, 15, 15, 15, 15, 15, 15

     0,   1.00e+00,   0.00e+00,   0.00e+00,   0.00e+00,   0.00e+00,   0.00e+00,   0.00e+00
   500,   8.54e-01,   4.57e-04,   1.44e-01,   1.55e-04,   6.47e-14,   1.23e-22,   6.44e-04
  1000,   7.30e-01,   3.90e-04,   2.65e-01,   2.89e-04,   2.53e-13,   2.28e-22,   2.44e-03
  1500,   6.23e-01,   3.33e-04,   3.66e-01,   4.02e-04,   2.98e-13,   3.18e-22,   5.20e-03
  2000,   5.32e-01,   2.85e-04,   4.49e-01,   5.00e-04,   3.30e-13,   3.95e-22,   8.77e-03
  2500,   4.55e-01,   2.43e-04,   5.18e-01,   5.83e-04,   3.55e-13,   4.61e-22,   1.30e-02
  3000,   3.88e-01,   2.08e-04,   5.75e-01,   6.54e-04,   3.74e-13,   5.17e-22,   1.78e-02
  3500,   3.32e-01,   1.77e-04,   6.21e-01,   7.14e-04,   3.88e-13,   5.65e-22,   2.30e-02
  4000,   2.83e-01,   1.52e-04,   6.59e-01,   7.66e-04,   4.00e-13,   6.06e-22,   2.86e-02
  4500,   2.42e-01,   1.29e-04,   6.88e-01,   8.10e-04,   4.09e-13,   6.41e-22,   3.45e-02
  5000,   2.07e-01,   1.11e-04,   7.11e-01,   8.48e-04,   4.15e-13,   6.71e-22,   4.06e-02
