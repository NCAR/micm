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
  integer, parameter :: dfactor = 1                               ! duplication factor to replicate the input data
  integer, parameter :: btime = 1, etime = 1, &                   ! begin & end index of each dimension of CAM FV output
                        dtime = etime-btime+1, &
                        blev = 1, elev = 32, &
                        dlev = elev-blev+1, &
                        blat = 96, elat = 97, &
                        dlat = elat-blat+1, &
                        blon = 130, elon = 134, &
                        dlon = elon-blon+1
  integer, parameter :: ntime = 1, nlev = 32, &                   ! For CAM FV 1-deg output
                        nlat = 192, nlon = 288
  integer, parameter :: kNumberOfGridCells = (elev-blev+1) * &    ! Number of grid cells with independent chemical reactions
                                             (elat-blat+1) * &
                                             (elon-blon+1) * &
                                             dfactor 
#else
  integer, parameter :: kNumberOfGridCells = 100                  ! Number of grid cells with independent chemical reactions
#endif
  integer, parameter :: masterproc = 0

end module constants
