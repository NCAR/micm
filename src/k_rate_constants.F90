module k_rateConst_module


use precision, only : r8

implicit none
private

public :: k_rateConst_register
public :: k_rateConst_init
public :: k_rateConst_run

! k_rateConst are computed at the beginning of the 
!   chemistry_box_solver time step.
!   They are not computed for internal steps of the
!   box-model time-step advancer
! k_rateConst will be thread-safe memory provided elsewhere.
! rate_constant_store will be an accessor to memory
! For now, it is allocated here. It is not thread safe

contains

  !---------------------------
  ! Register the number of k_rateConst values for the run
  !---------------------------
  subroutine k_rateConst_register(nkReact)
      
    integer, parameter :: nkReact_set = 3
    integer :: nkReact

    nkReact=nkReact_set

  end  subroutine k_rateConst_register

  !---------------------------
  ! Register the number of k_rateConst values for the run
  !---------------------------
  subroutine k_rateConst_init(nkReact, k_rateConst)
    integer, intent(in) :: nkReact 
    real(r8), pointer, intent(inout) :: k_rateConst(:)
  end  subroutine k_rateConst_init

  !---------------------------
  ! Compute k_rateConst, given M, P, T
  ! Execute once for the chemistry-time-step advance
  ! Not called from the solver
  !---------------------------
  subroutine k_rateConst_run(k_rateConst)
  
    use external_fields, M => mass_density, P => pressure, T=>temperature
  real(r8) :: t_inverse
  real(r8),pointer :: k_rateConst(:)
  
  t_inverse = 1/T

!>>FromCafe
! Rate Constants
! Y0_a
k_rateConst(1) = 0.04
! Y1_Y2_M_b
k_rateConst(2) = 1e+4
! Y1_Y1_a
k_rateConst(3) = 1.5e7 * exp(0 * t_inverse)
!<<FromCafe

  print*,'rate constants'
  print*,k_rateConst
  
  end subroutine k_rateConst_run
  
  
end module k_rateConst_module
