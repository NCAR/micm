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

    .. tab-item:: OpenAtmos Configuration reading

        .. raw:: html

            <div class="download-div">
              <a href="../_static/tutorials/rate_constants_user_defined.zip" download>
                <button class="download-button">Download zip configuration</button>
              </a>
            </div>

        .. literalinclude:: ../../../test/tutorial/test_rate_constants_user_defined_with_config.cpp
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

             Process r7 = Process::create()
                             .reactants({ f })
                             .products({ yields(g, 1) })
                             .rate_constant(TunnelingRateConstant({ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 }))
                             .phase(gas_phase);

            + Process r8 = Process::create()
            +                 .reactants({ c })
            +                 .products({ yields(g, 1) })
            +                 .rate_constant(UserDefinedRateConstant({.label_="my rate"}))
            +                 .phase(gas_phase);

            + Process r9 = Process::create()
            +                 .products({ yields(a, 1) })
            +                 .rate_constant(UserDefinedRateConstant({.label_="my emission rate"}))
            +                 .phase(gas_phase);

            + Process r10 = Process::create()
            +                 .reactants({ b })
            +                 .rate_constant(UserDefinedRateConstant({.label_="my loss rate"}))
            +                 .phase(gas_phase);

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

          +double photo_rate = 1e-10;
          +double emission_rate = 1e-20;
          +double loss = emission_rate * 1e-3;
          +// these rates are constant through the simulation
          +state.SetCustomRateParameter("my emission rate", emission_rate);
          +state.SetCustomRateParameter("my loss rate", loss);
           // solve for ten iterations
           for (int i = 0; i < 10; ++i)
           {
             // Depending on how stiff the system is
             // the solver integration step may not be able to solve for the full time step
             // so we need to track how much time the solver was able to integrate for and continue
             // solving until we finish
             double elapsed_solve_time = 0;
          +  state.SetCustomRateParameter("my photolysis rate", photo_rate);

             while (elapsed_solve_time < time_step)
             {
               auto result = solver.Solve(time_step - elapsed_solve_time, state);
               elapsed_solve_time = result.final_time_;
               // std::cout << "solver state: " << StateToString(result.state_) << std::endl;
               state.variables_[0] = result.result_.AsVector();
             }

             print_state(time_step * (i + 1), state);
          +   photo_rate *= 1.5;
           }

    .. tab-item:: OpenAtmos Configuration reading

        In this case, you only need to add the configuration to the reactions.json file in the configuration directory.
        When reading in from a configuration file, the loss, emissions, and photolysis rates are prefixed with
        ``LOSS.``, ``EMIS.``, and ``PHOTO.``. This differs slightly from defining the API by hand.


        .. code-block:: diff

          +double photo_rate = 1e-10;
          +double emission_rate = 1e-20;
          +double loss = emission_rate * 1e-3;
          +// these rates are constant through the simulation
          +state.SetCustomRateParameter("EMIS.my emission rate", emission_rate);
          +state.SetCustomRateParameter("LOSS.my loss rate", loss);
           // solve for ten iterations
           for (int i = 0; i < 10; ++i)
           {
             // Depending on how stiff the system is
             // the solver integration step may not be able to solve for the full time step
             // so we need to track how much time the solver was able to integrate for and continue
             // solving until we finish
             double elapsed_solve_time = 0;
             +state.SetCustomRateParameter("PHOTO.my photolysis rate", photo_rate);

             while (elapsed_solve_time < time_step)
             {
               auto result = solver.Solve(time_step - elapsed_solve_time, state);
               elapsed_solve_time = result.final_time_;
               // std::cout << "solver state: " << StateToString(result.state_) << std::endl;
               state.variables_[0] = result.result_.AsVector();
             }

             print_state(time_step * (i + 1), state);
             +photo_rate *= 1.5;
           }

And this is final output. Notice that the concentration of G ends up much higher than in 
the :ref:`Rate constants` tutorial's result.

.. csv-table:: Table Title
   :header: "time", "A", "B", "C", "D", "E", "F", "G"
   :widths: 10, 15, 15, 15, 15, 15, 15, 15

   "0", "1.00e+00", "0.00e+00", "0.00e+00", "0.00e+00", "0.00e+00", "0.00e+00", "0.00e+00"
   "500", "3.18e-09", "3.66e-09", "9.83e-01", "3.88e-14", "1.41e-03", "2.02e-13", "7.92e-03"
   "1000", "1.14e-14", "1.31e-14", "9.66e-01", "1.39e-19", "1.40e-03", "7.24e-19", "1.64e-02"
   "1500", "7.27e-20", "6.40e-20", "9.49e-01", "6.53e-25", "1.39e-03", "3.19e-24", "2.48e-02"
   "2000", "3.17e-20", "1.70e-20", "9.33e-01", "1.55e-25", "1.38e-03", "5.92e-25", "3.30e-02"
   "2500", "3.17e-20", "1.70e-20", "9.17e-01", "1.55e-25", "1.37e-03", "5.92e-25", "4.11e-02"
   "3000", "3.17e-20", "1.70e-20", "9.01e-01", "1.55e-25", "1.36e-03", "5.92e-25", "4.90e-02"
   "3500", "3.17e-20", "1.70e-20", "8.85e-01", "1.55e-25", "1.34e-03", "5.92e-25", "5.68e-02"
   "4000", "3.17e-20", "1.70e-20", "8.70e-01", "1.55e-25", "1.33e-03", "5.92e-25", "6.44e-02"
   "4500", "3.17e-20", "1.70e-20", "8.55e-01", "1.55e-25", "1.32e-03", "5.92e-25", "7.20e-02"
   "5000", "3.17e-20", "1.70e-20", "8.40e-01", "1.55e-25", "1.31e-03", "5.92e-25", "7.94e-02"

