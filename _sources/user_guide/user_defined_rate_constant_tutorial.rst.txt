.. _User defined rate constants:

User-Defined Rate Constants
###########################

This tutorial extends the :ref:`Rate constants` tutorial. That one showed how to use every type of rate constant
except the user-defined ones. The difference is that a user-defined rate constant must be updated by the user,
whereas the other rate constants update in response to a changing state. 

Internal to MICM, user-defined rate constants provide the ability to represent 
processes like emissions, loss, and photolysis that have rate constants that are calculated externally from MICM. These can be static and the same for each time step, or dynamically updated 
during the simulation.

To show how this works, we'll add one photolysis reaction, one first-order loss reaction, and one emission

.. math::

  A &\longrightarrow B, &k_{1, \mathrm{arrhenius}} \\
  B &\longrightarrow C (\mathrm{alkoxy\ products}) + D (\mathrm{nitrate\ products}), &k_{2, \mathrm{branched}} \\
  C &\longrightarrow E, &k_{3, \mathrm{surface}} \\
  D &\longrightarrow 2F, &k_{4, \mathrm{ternary\ chemical\ activation}} \\
  2E &\longrightarrow G, &k_{5, \mathrm{troe}} \\
  F &\longrightarrow G, &k_{6, \mathrm{tunneling}} \\
  C &\longrightarrow G, &k_{7, \mathrm{photolysis}} \\
  &\longrightarrow A, &k_{8, \mathrm{emission}} \\
  B &\longrightarrow, &k_{9, \mathrm{loss}} \\


If you're looking for a copy and paste, choose
the appropriate tab below and be on your way! Otherwise, stick around for a line by line explanation.

.. tab-set::

  .. tab-item:: Build the Mechanism with the API

    .. literalinclude:: ../../../test/tutorial/test_rate_constants_user_defined_by_hand.cpp
      :language: cpp

Line-by-line explanation
------------------------

Adding the custom rate constant is quite simple. Include the header file:

.. code-block:: diff

  #include <micm/process/arrhenius_rate_constant.hpp>
  #include <micm/process/branched_rate_constant.hpp>
  #include <micm/process/surface_rate_constant.hpp>
  #include <micm/process/ternary_chemical_activation_rate_constant.hpp>
  #include <micm/process/troe_rate_constant.hpp>
  #include <micm/process/tunneling_rate_constant.hpp>
  + #include <micm/process/user_defined_rate_constant.hpp>
  #include <micm/solver/rosenbrock.hpp>


Then setup the reaction which will use this rate constant:

.. tab-set::

  .. tab-item:: Build the Mechanism with the API

    .. code-block:: diff

      Process r7 = Process::Create()
                      .SetReactants({ f })
                      .SetProducts({ Yields(g, 1) })
                      .SetRateConstant(TunnelingRateConstant({ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 }))
                      .SetPhase(gas_phase);

      + Process r8 = Process::Create()
      +                 .SetReactants({ c })
      +                 .SetProducts({ Yields(g, 1) })
      +                 .SetRateConstant(UserDefinedRateConstant({.label_="my rate"}))
      +                 .SetPhase(gas_phase);

      + Process r9 = Process::Create()
      +                 .SetProducts({ Yields(a, 1) })
      +                 .SetRateConstant(UserDefinedRateConstant({.label_="my emission rate"}))
      +                 .SetPhase(gas_phase);

      + Process r10 = Process::Create()
      +                 .SetReactants({ b })
      +                 .SetRateConstant(UserDefinedRateConstant({.label_="my loss rate"}))
      +                 .SetPhase(gas_phase);

      auto chemical_system = System(micm::SystemParameters{ .gas_phase_ = gas_phase });
      - auto reactions = std::vector<micm::Process>{ r1, r2, r3, r4, r5, r6, r7 };
      + auto reactions = std::vector<micm::Process>{ r1, r2, r3, r4, r5, r6, r7, r8, r9, r10 };


  .. tab-item:: OpenAtmos Configuration reading

    In this case, you only need to add the configuration to the reactions.json file in the configuration directory.

    .. code-block:: diff

      + {
      +   "type": "PHOTOLYSIS",
      +   "reactants": {
      +     "C": {}
      +   },
      +   "products": {
      +     "G": {}
      +   },
      +   "MUSICA name": "my photolysis rate"
      + },
      + {
      +   "type": "FIRST_ORDER_LOSS",
      +   "species": "B",
      +   "MUSICA name": "my loss rate"
      + },
      + {
      +   "type": "EMISSION",
      +   "species": "A",
      +   "MUSICA name": "my emission rate"
      + }


Finally, set and upate the rate constants as needed:


.. tab-set::

  .. tab-item:: Build the Mechanism with the API

    .. code-block:: diff

      + double photo_rate = 1e-10;
      + double emission_rate = 1e-20;
      + double loss = emission_rate * 1e-3;
      + // these rates are constant through the simulation
      + state.SetCustomRateParameter("my emission rate", emission_rate);
      + state.SetCustomRateParameter("my loss rate", loss);
        // solve for ten iterations
        for (int i = 0; i < 10; ++i)
        {
          // Depending on how stiff the system is
          // the solver integration step may not be able to solve for the full time step
          // so we need to track how much time the solver was able to integrate for and continue
          // solving until we finish
          double elapsed_solve_time = 0;
      +   state.SetCustomRateParameter("my photolysis rate", photo_rate);
          solver.CalculateRateConstants(state);

          while (elapsed_solve_time < time_step)
          {
            auto result = solver.Solve(time_step - elapsed_solve_time, state);
            elapsed_solve_time = result.final_time_;
          }

          print_state(time_step * (i + 1), state);
      +   photo_rate *= 1.5;
        }

  .. tab-item:: OpenAtmos Configuration reading

    In this case, you only need to add the configuration to the reactions.json file in the configuration directory.
    When reading in from a configuration file, the loss, emissions, and photolysis rates are prefixed with
    ``LOSS.``, ``EMIS.``, and ``PHOTO.``. This differs slightly from defining the API by hand.

    .. code-block:: diff

      + double photo_rate = 1e-10;
      + double emission_rate = 1e-20;
      + double loss = emission_rate * 1e-3;
      + // these rates are constant through the simulation
      + state.SetCustomRateParameter("EMIS.my emission rate", emission_rate);
      + state.SetCustomRateParameter("LOSS.my loss rate", loss);

        // solve for ten iterations
        for (int i = 0; i < 10; ++i)
        {
          // Depending on how stiff the system is
          // the solver integration step may not be able to solve for the full time step
          // so we need to track how much time the solver was able to integrate for and continue
          // solving until we finish
          double elapsed_solve_time = 0;
      +   state.SetCustomRateParameter("PHOTO.my photolysis rate", photo_rate);
          solver.CalculateRateConstants(state);

          while (elapsed_solve_time < time_step)
          {
            auto result = solver.Solve(time_step - elapsed_solve_time, state);
            elapsed_solve_time = result.final_time_;
          }

          print_state(time_step * (i + 1), state);
      +   photo_rate *= 1.5;
        }

And this is final output. Notice that the concentration of G ends up much higher than in 
the :ref:`Rate constants` tutorial's result.

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
