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

We'll setup and solve a fake chemical system with 4 species and 6 reactions, 

.. math::

  A &\longrightarrow B, &k_{1, \mathrm{arrhenius}} \\
  B &\longrightarrow C (\mathrm{alkoxy\ products}) + D (\mathrm{nitrate\ products}), &k_{2, \mathrm{branched}} \\
  2C &\longrightarrow A, &k_{3, \mathrm{surface}} \\
  D &\longrightarrow 2B, &k_{4, \mathrm{ternary\ chemical\ activation}} \\
  A &\longrightarrow C, &k_{5, \mathrm{troe}} \\
  B &\longrightarrow D, &k_{6, \mathrm{tunneling}} \\


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
    :lines: 1-11

After that, we'll use the ``micm`` namespace and setup a template alias so that we can instantiate the 
rosenbrock solver.

  .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
    :language: cpp
    :lines: 13-20

Next we have to define the chemical system we want to solve. 

.. tabs::

    .. tab:: Build the Mechanism with the API

        To do this by hand, we have to define all of the chemical species in the system. This allows us to set
        any properties of the species that may be necessary for rate constanta calculations, like molecular weights 
        and diffusion coefficients for the surface reaction. 

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_by_hand.cpp
          :language: cpp
          :lines: 25-31

        Products are created with a :cpp:type:`micm::Yield` :cpp:func:`micm::yields`

    .. tab:: OpenAtmos Configuration reading

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_no_user_defined_with_config.cpp
          :language: cpp