! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Stub functions for the mechanism rate constants

!> Mechanism-specific rate constant functions
module rate_constants_utility

  use micm_environment,                only : environment_t
  use micm_rate_constant_arrhenius,    only : rate_constant_arrhenius_t
  use micm_rate_constant_photolysis,   only : rate_constant_photolysis_t
  use micm_rate_constant_ternary_chemical_activation,                         &
      only : rate_constant_ternary_chemical_activation_t
  use micm_rate_constant_troe,                                                &
      only : rate_constant_troe_t
  use micm_rate_constant_wennberg_alkoxy,                                     &
      only : rate_constant_wennberg_alkoxy_t
  use micm_rate_constant_wennberg_nitrate,                                    &
      only : rate_constant_wennberg_nitrate_t
  use micm_rate_constant_wennberg_tunneling,                                  &
      only : rate_constant_wennberg_tunneling_t
  use musica_constants,                only : musica_dk

  implicit none
  private

  public :: calculate_rate_constants

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the rate constant for each reaction
  subroutine calculate_rate_constants( rate_constants, environment )

    !> Rate constant for each reaction [(molec cm-3)^(n-1) s-1]
    real(kind=musica_dk), intent(out) :: rate_constants(:)
    !> Environmental state
    type(environment_t),  intent(in)  :: environment

    type( rate_constant_arrhenius_t                   ) :: arrhenius
    type( rate_constant_photolysis_t                  ) :: photolysis
    type( rate_constant_ternary_chemical_activation_t ) :: ternary_chemical_activation
    type( rate_constant_troe_t                        ) :: troe
    type( rate_constant_wennberg_alkoxy_t             ) :: wennberg_alkoxy
    type( rate_constant_wennberg_nitrate_t            ) :: wennberg_nitrate
    type( rate_constant_wennberg_tunneling_t          ) :: wennberg_tunneling

  end subroutine calculate_rate_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module rate_constants_utility
