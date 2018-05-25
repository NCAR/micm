program micm_driver

!--------------------------------------------------------------------------------
! Program which prototypes the MICM driver
!--------------------------------------------------------------------------------
use external_fields,     only: set_externals
use chemSolve,          only: chemSolve_register,  chemSolve_init, chemSolve_run
use k_rateConst_module, only: k_rateConst_register, k_rateConst_init, k_rateConst_run

use precision, only : r8

implicit none

integer :: nSpecies                    ! Number of chemical species in run
integer :: nkReact                     ! Number of k reactions in run
real(r8), pointer :: k_rateConst(:)

real(r8),pointer :: vmr_init(:)
real(r8),pointer :: vmr_curr(:)

! this system is unstable, so 3e+1 fails
real(r8) :: time_step_size = 2e+1_r8 ! seconds

! convergence criteria will have to be set somewhere and passed to ode solver.
real(r8),pointer :: RelTol(:)
real(r8),pointer :: AbsTol(:)
integer :: ode_retcode

!-----------------------------------------------
! Register the chemistry packages
!-----------------------------------------------

call chemSolve_register(nSpecies)
call k_rateConst_register(nkReact)


!-----------------------------------------------
! Allocate the local variables  (This will be done via CPF?)
!-----------------------------------------------
allocate (vmr_init(nSpecies))
allocate (vmr_curr(nSpecies))
allocate (AbsTol(nSpecies))
allocate (RelTol(nSpecies))
allocate(k_rateConst(nkReact))

!-----------------------------------------------
! Initialize the chemistry packages
!-----------------------------------------------

call k_rateConst_init(nkReact, k_rateConst)
call chemSolve_init(absTol, relTol)

!-----------------------------------------------
! Explicitly specify the data which will come from the  cpf
!-----------------------------------------------

  vmr_init(:) =(/1._r8, 0._r8, 0._r8/)
  print *, 'initial value'
  print *, vmr_init


!-----------------------------------------------
! Simulate the XML file which CCPP will use to drive the model
!-----------------------------------------------
! Only called at beginnning
  call k_rateConst_run(k_rateConst)

  call  chemSolve_run(nkReact, vmr_init, time_step_size, k_rateConst, AbsTol, RelTol, vmr_curr, ode_retcode)

! data back to cpf

  print *, 'final_value'
  print *, vmr_curr

!-----------------------------------
! deallocate variables
!-----------------------------------

  deallocate (vmr_init, vmr_curr, k_rateConst)

end program micm_driver
