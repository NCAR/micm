module forcing_and_jacobian

use precision, only : r8
use chemistry_specification, only: n_reactions

implicit none

contains

  !---------------------------
  ! Compute reaction rates, given vmr of each species and rate constants that have been computed elsewhere
  !---------------------------
  function compute_rates(vmr)
  
    use external_fields, M =>mass_density
    use rate_constants, only: get_rate_constants
  
    real(r8), intent(in) ::  vmr(:)           ! volume mixing ratios of each component in order
    real(r8) ::  rateConstants(n_reactions)  ! rate_constants for each reaction
    real(r8) ::  rates(n_reactions)  ! rate_constants for each reaction
    real(r8) ::  compute_rates(n_reactions)   ! rates for each reaction (sometimes called velocity of reaction)
  
    ! this could be access to memory that has been allocated by cpf, or host model.
    rateConstants = get_rate_constants()
    
   include 'rates'

    compute_rates = rates
  
  end function compute_rates
  
  
  !---------------------------
  ! Compute time rate of change of each molecule (vmr) given reaction rates
  !---------------------------
  function force(vmr)
  
    real(r8), intent(in)::  vmr(:)   ! volume mixing ratios of each component in order
    real(r8) ::  force(size(vmr))    ! rate of change of each molecule
    real(r8) ::  rates(n_reactions)  ! rates of each reaction
  
    rates = compute_rates(vmr)
 
   include 'forcing'
 
  end function force
  
  
  !---------------------------
  ! Compute sensitivity of molecular forcing to each vmr (derivative of force w.r.t. each vmr)
  !---------------------------
  function jac(vmr)
  
    use rate_constants, only: get_rate_constants
    use external_fields, M=>  mass_density

    real(r8), intent(in)::  vmr(:)              ! volume mixing ratios of each component in order
    real(r8) :: jac(size(vmr),size(vmr))   ! sensitivity of forcing to changes in each vmr
    real(r8) ::  rateConstants(n_reactions)    ! forcing (rate of change) of each component in order
  
    jac(:,:) = 0.
  
    rateConstants = get_rate_constants()

    include 'jacobian'

  end function jac

end module forcing_and_jacobian
