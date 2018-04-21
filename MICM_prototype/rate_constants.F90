module rate_constants


use precision, only : r8
use chemistry_specification, only: n_reactions

implicit none

! rate_constants are computed at the beginning of the 
!   chemistry_box_solver time step.
!   They are not computed for internal steps of the
!   box-model time-step advancer
! rate_constants will be thread-safe memory provided elsewhere.
! rate_constant_store will be an accessor to memory
! For now, it is allocated here. It is not thread safe

real(r8) :: rateConstants(n_reactions)

contains

  !---------------------------
  ! Compute rate_constants, given M, P, T
  ! Execute once for the chemistry-time-step advance
  ! Not called from the solver
  !---------------------------
  subroutine compute_rate_constants()
  
    use external_fields, M => mass_density, P => pressure, T=>temperature
    !use rate_constant_store, only: get_stored_rate_constants
  real(r8) :: t_inverse
  
  t_inverse = 1/T
  
  include 'rateconstants'

  print*,'rate constants'
  print*,rateConstants
  
  end subroutine compute_rate_constants
  
  
  !---------------------------
  ! Provide rate_constants from store
  ! This may be called multiple times from the solver
  !---------------------------
  function get_rate_constants()
  real(r8), dimension(n_reactions) :: get_rate_constants
  
    ! rate_constants memory is allocated elsewhere
    get_rate_constants = rateConstants
 
  end function get_rate_constants
  
  
end module rate_constants
