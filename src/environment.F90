! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \todo rewrite micm_environment module to follow musica style conventions
module micm_environment

  use kinetics_utilities,              only : number_of_photolysis_reactions
  use musica_constants,                only : musica_dk
  use constants,                       only : ncell=>kNumberOfGridCells

  type environment_t
    !> Number density of air [# cm-3]
    real(musica_dk) :: number_density_air(ncell)
    !> Temperature [K]
    real(musica_dk) :: temperature(ncell)
    !> Pressure [Pa]
    real(musica_dk) :: pressure(ncell)
    !> Aerosol surface area density \bug what are the units for surface area density?
    real(musica_dk) :: aerosol_surface_area_density(ncell,4)
    !> Aerosol diameter \bug what are the units for aerosol diameter?
    real(musica_dk) ::  aerosol_diameter(ncell,4)
    !> Water vapor volume mixing ratio [mol mol-1]
    real(musica_dk) :: h2ovmr(ncell)
    !> Molecular oxygen number density [# cm-3]
    real(musica_dk) :: o2_number_density(ncell)
    !> Photolysis rate constants [s-1]
    real(musica_dk) :: photolysis_rate_constants(ncell,number_of_photolysis_reactions)
  end type environment_t

end module micm_environment
