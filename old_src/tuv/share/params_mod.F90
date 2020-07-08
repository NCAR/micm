module params_mod

  use phot_kind_mod, only: rk => kind_phot

  implicit none

  public

  real(rk), parameter :: m2km = .001_rk            ! meters to km
  real(rk), parameter :: ppm2vmr = 1.e-6_rk        ! ppm to vmr
  real(rk), parameter :: km2cm = 1.e5_rk           ! km to centimeters
  real(rk), parameter :: m2s   = 60._rk            ! minutes to seconds

  real(rk), parameter :: pi = 3.1415926535898_rk
  real(rk), parameter :: radius = 6.371E+3_rk      ! km
  real(rk), parameter :: hc = 6.626068E-34_rk * 2.99792458E8_rk ! (J sec)*(m/sec) ==> watts sec m
  real(rk), parameter :: largest=1.E+36_rk
  real(rk), parameter :: kboltz= 1.38064852e-16_rk ! boltzmann constant (erg/K)
  real(rk), parameter :: R=2.8704e6_rk             ! gas constant (erg/g/K)
  real(rk), parameter :: g=980.616_rk              ! grav acceleration (cm/sec2)


  real(rk), parameter :: pzero = +10._rk/largest
  real(rk), parameter :: nzero = -10._rk/largest

  real(rk), parameter :: precis = 1.e-7_rk

  character(len=256) :: input_data_root = 'NOT_SET'

  ! input unit number
  integer, parameter :: kin=12

  ! delta for adding points at beginning or end of data grids
  real(rk), parameter :: deltax = 1.E-5_rk

contains

  ! quiet NaN
 function qnan() result(x)
    USE,INTRINSIC :: IEEE_ARITHMETIC
    real(rk) :: x
    x = IEEE_VALUE(x,IEEE_QUIET_NAN)
  end function qnan
 
end module params_mod
