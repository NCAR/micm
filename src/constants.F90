! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> Physical constants
module constants

  use musica_constants, only : dk => musica_dk

  real(kind=dk), parameter :: kBoltzmann = 1.38065e-23_dk         ! J/K/molecule
  real(kind=dk), parameter :: kAvagadro  = 6.02214076e23_dk       ! molecules / mole
  integer, parameter :: VLEN = 128                                ! vector length for GPU kernels
  integer, parameter :: STREAM0 = 0                               ! stream ID for async GPU kernels
#ifdef USE_NETCDF
  integer, parameter :: ntime = 1, nlev = 32, &
                        nlat = 192, nlon = 288                    ! For CAM FV 1-deg output
  integer, parameter :: kNumberOfGridCells = nlev*nlat*nlon       ! Number of grid cells with independent chemical reactions
#else
  integer, parameter :: kNumberOfGridCells = 100                  ! Number of grid cells with independent chemical reactions
#endif
!  integer, parameter :: masterproc = 0
!  integer :: beg_grid, end_grid

end module constants
