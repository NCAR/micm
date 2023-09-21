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
in the OpenAtmos format. In this tutorial, we will do both. If you're looking for a copy and paste, choose
the appropriate tab below and be on your way! Otherwise, stick around for a line by line explanation.


.. tabs::

    .. tab:: Build the Mechanism with the API

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp

    .. tab:: OpenAtmos Configuration reading

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

After that, we'll use the ``micm`` namespace and setup a template alias so that we can instantiate the 
rosenbrock solver.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 15-22

To create a :cpp:class:`micm::RosenbrockSolver`, we have to define a chemical system (:cpp:class:`micm::System`)
and our reactions, which will be a vector of :cpp:class:`micm::Process` We will use the species to define these.

.. tabs::

    .. tab:: Build the Mechanism with the API

        To do this by hand, we have to define all of the chemical species in the system. This allows us to set
        any properties of the species that may be necessary for rate constanta calculations, like molecular weights 
        and diffusion coefficients for the surface reaction.  We will also put these species into the gas phase.

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 56-67

        Now that we have a gas phase and our species, we can start building the reactions. Two things to note are that
        stoichiemtric coefficients for reactants are represented by repeating that product as many times as you need.
        To specify the yield of a product, we've created a typedef :cpp:type:`micm::Yield` 
        and a function :cpp:func:`micm::yields` that produces these.

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 69-133
        
        And finally we define our chemical system and reactions

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 135-136

    .. tab:: OpenAtmos Configuration reading

        After defining a valid OpenAtmos configuration with reactions that ``micm`` supports, configuring the chemical
        system and the processes is as simple as using the :cpp:class:`micm::SolverConfig` class

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_with_config.cpp
          :language: cpp
          :lines: 57-69

Now that we have a chemical system and a list of reactions, we can create the RosenbrockSolver.
There are several ways to configure the solver. Here we are using a three stage solver. More options
can be found in the :cpp:class:`micm::RosenbrockSolverParameters`

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 138-140

The rosenbrock solver will provide us a state, which we can use to set the concentrations,
custom rate parameters, and temperature and pressure

.. tabs::

    .. tab:: Build the Mechanism with the API

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 141-155

    .. tab:: OpenAtmos Configuration reading

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_with_config.cpp
          :language: cpp
          :lines: 75-93


Finally, we are ready to pick a timestep ans solve the system.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 157-183


This is the output:


+-------+------------+------------+------------+------------+------------+------------+------------+
| time  |     A      |     B      |     C      |     D      |     E      |     F      |     G      |
+=======+============+============+============+============+============+============+============+
|   0   | 1.00e+00   | 0.00e+00   | 0.00e+00   | 0.00e+00   | 0.00e+00   | 0.00e+00   | 0.00e+00   |
+-------+------------+------------+------------+------------+------------+------------+------------+
|  500  | 3.22e-09   | 3.70e-09   | 9.67e-01   | 3.92e-14   | 1.38e-03   | 2.04e-13   | 7.69e-03   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 1000  | 1.15e-14   | 1.33e-14   | 9.35e-01   | 1.40e-19   | 1.34e-03   | 7.31e-19   | 1.56e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 1500  | 4.14e-20   | 4.76e-20   | 9.06e-01   | 5.04e-25   | 1.29e-03   | 2.62e-24   | 2.30e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 2000  | 1.48e-25   | 1.71e-25   | 8.78e-01   | 1.81e-30   | 1.26e-03   | 9.40e-30   | 3.00e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 2500  | 5.32e-31   | 6.12e-31   | 8.52e-01   | 6.47e-36   | 1.22e-03   | 3.37e-35   | 3.65e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 3000  | 1.91e-36   | 2.19e-36   | 8.27e-01   | 2.32e-41   | 1.18e-03   | 1.21e-40   | 4.27e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 3500  | 6.84e-42   | 7.86e-42   | 8.04e-01   | 8.32e-47   | 1.15e-03   | 4.33e-46   | 4.85e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 4000  | 2.45e-47   | 2.82e-47   | 7.82e-01   | 2.98e-52   | 1.12e-03   | 1.55e-51   | 5.40e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 4500  | 8.80e-53   | 1.01e-52   | 7.61e-01   | 1.07e-57   | 1.09e-03   | 5.57e-57   | 5.92e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+
| 5000  | 3.16e-58   | 3.63e-58   | 7.42e-01   | 3.84e-63   | 1.06e-03   | 2.00e-62   | 6.41e-02   |
+-------+------------+------------+------------+------------+------------+------------+------------+