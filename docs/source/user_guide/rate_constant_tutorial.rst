.. _Rate constants:

Rate Constants (except user defined ones)
#########################################

MICM supports a subset of the rate constants defined as part of the 
`OpenAtmos Mechanism Configuration <https://open-atmos.github.io/MechanismConfiguration/reactions/index.html>`_
We will be adding more in the future. At present, supported rate constants are:

- :cpp:class:`micm::ArrheniusRateConstant`
- :cpp:class:`micm::BranchedRateConstant`
- :cpp:class:`micm::SurfaceRateConstant`
- :cpp:class:`micm::TernaryChemicalActivationRateConstant`
- :cpp:class:`micm::TroeRateConstant`
- :cpp:class:`micm::TunnelingRateConstant`
- :cpp:class:`micm::UserDefinedRateConstant`

This tutorial covers all but the last one. See the :ref:`User defined rate constants` tutorial for examples and use
cases on that.
