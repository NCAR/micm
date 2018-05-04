module k_rate_const_module


use precision, only : r8

implicit none
private

public :: k_rate_const_register
public :: k_rate_const_init
public :: k_rate_const_run

! k_rate_const are computed at the beginning of the 
!   chemistry_box_solver time step.
!   They are not computed for internal steps of the
!   box-model time-step advancer
! k_rate_const will be thread-safe memory provided elsewhere.
! rate_constant_store will be an accessor to memory
! For now, it is allocated here. It is not thread safe

contains

  !---------------------------
  ! Register the number of k_rate_const values for the run
  !---------------------------
  subroutine k_rate_const_register(nkReact)
      
    integer, parameter :: nkReact_set = 3
    integer :: nkReact

    nkReact=nkReact_set

  end  subroutine k_rate_const_register

  !---------------------------
  ! Register the number of k_rate_const values for the run
  !---------------------------
  subroutine k_rate_const_init(nkReact, k_rate_const)
    integer, intent(in) :: nkReact 
    real(r8), pointer, intent(inout) :: k_rate_const(:)
  end  subroutine k_rate_const_init

  !---------------------------
  ! Compute k_rate_const, given M, P, T
  ! Execute once for the chemistry-time-step advance
  ! Not called from the solver
  !---------------------------
  subroutine k_rate_const_run(k_rate_const)
  
    use external_fields, M => mass_density, P => pressure, T=>temperature
    !use rate_constant_store, only: get_stored_k_rate_const
  real(r8) :: t_inverse
  real(r8),pointer :: k_rate_const(:)
  
  t_inverse = 1/T
  ! *** CALL OR USE ASSOC THIS INCLUDE
  include 'rateconstants'

  print*,'rate constants'
  print*,k_rate_const
  
  end subroutine k_rate_const_run
  
  
end module k_rate_const_module
