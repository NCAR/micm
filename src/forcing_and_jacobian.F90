module forcing_and_jacobian

use precision, only : r8

implicit none

contains

  !---------------------------
  ! Compute reaction rates, given vmr of each species and rate constants that have been computed elsewhere
  !---------------------------
  subroutine compute_rates(nkReact, vmr, k_rate_const,  k_rates)
  
    use external_fields, M =>mass_density
  
    integer, intent(in) :: nkReact
    real(r8), intent(in) ::  vmr(:)           ! volume mixing ratios of each component in order
    real(r8), intent(in) ::  k_rate_const(nkReact)   ! rates constants for each reaction
    real(r8), intent(out) ::  k_rates(nkReact)   ! rates for each reaction (sometimes called velocity of reaction)
  
! Rates
! Y0_a
k_rates(1) = k_rate_const(1) * vmr(1)
! Y1_Y2_M_b
k_rates(2) = k_rate_const(2) * vmr(2) * vmr(3) * M
! Y1_Y1_a
k_rates(3) = k_rate_const(3) * vmr(2) * vmr(2)

  end subroutine compute_rates
  
  
  !---------------------------
  ! Compute time rate of change of each molecule (vmr) given reaction rates
  !---------------------------
  subroutine calc_force(nkReact, vmr, k_rate_const, force)
  
    integer, intent(in)::  nkReact
    real(r8), intent(in)::  vmr(:)   ! volume mixing ratios of each component in order
    real(r8), intent(in)::  k_rate_const(:)   ! rates constants for each reaction
    real(r8), intent(out) ::  force(size(vmr))    ! rate of change of each molecule

    real(r8) ::  k_rates(nkReact)  ! rates of each reaction
  
    call compute_rates(nkReact,vmr,k_rate_const, k_rates)
 
! testing
 force(1) = (-1) * k_rates(1) + (1) * k_rates(2)
 force(2) = (1) * k_rates(1) + (-1) * k_rates(2) + (-1) * k_rates(3) + (-1) * k_rates(3)
 force(3) = (0) * k_rates(2) + (2) * k_rates(3)

  end subroutine calc_force
  
  
  !---------------------------
  ! Compute sensitivity of molecular forcing to each vmr (derivative of force w.r.t. each vmr)
  !---------------------------
  function jac(vmr, k_rate_const)
  
    use external_fields, M=>  mass_density

    real(r8), intent(in)::  vmr(:)              ! volume mixing ratios of each component in order
    real(r8) :: jac(size(vmr),size(vmr))   ! sensitivity of forcing to changes in each vmr
    real(r8),pointer ::  k_rate_const(:)    ! forcing (rate of change) of each component in order
  
    jac(:,:) = 0.
  
    include 'jacobian'

  end function jac

end module forcing_and_jacobian
