! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> Physical constants
module constants

  real, parameter :: kBoltzmann = 1.38065e-23                     ! J/K/molecule
  real, parameter :: kAvagadro  =  6.02214076e23                  ! molecules / mole
  integer, parameter :: VLEN = 128                                ! vector length for GPU kernels
  integer, parameter :: STREAM0 = 0                               ! stream ID for async GPU kernels
#ifdef USE_NETCDF
  integer, parameter :: ntime = 1, nlev = 32, &
                        nlat = 192, nlon = 288                    ! For CAM FV 1-deg output
  integer, parameter :: kNumberOfGridCells = nlev*nlat*nlon       ! Number of grid cells with independent chemical reactions
  !> For Chapman mechanism
  character(len=128), parameter :: cam_photo_reaction_names = &   ! Corresponding MICM photolysis reaction names in CAM output
                                   (/'jo2_a','jo3_a','jo3_b'/)    
#else
  integer, parameter :: kNumberOfGridCells = 100                  ! Number of grid cells with independent chemical reactions
#endif

end module constants
