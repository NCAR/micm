program host_model_simulator

!--------------------------------------------------------------------------------
! Program which simulates the host model for MICM
!--------------------------------------------------------------------------------
use chemSolve,          only: chemSolve_register,  chemSolve_init, chemSolve_run
use k_rateConst_module, only: k_rateConst_register, k_rateConst_init, k_rateConst_run
use precision,          only: r8

implicit none

integer           :: nSpecies               ! Number of chemical species in run
integer           :: nkReact                ! Number of k reactions in run
real(r8), pointer :: k_rateConst(:)         ! K rate constants

real(r8),pointer  :: vmr_init(:)            ! Initial VMR
real(r8),pointer  :: vmr_curr(:)            ! Current VMR

real(r8)          :: timeStepSize = 2e+1_r8 ! seconds - this system is unstable, so 3e+1 fails

! convergence criteria will have to be set somewhere(Cafe) and passed to ode solver.
real(r8),pointer  :: relTol(:)              ! Relative tolerance
real(r8),pointer  :: absTol(:)              ! Absolute tolerance

integer           :: ierr                   ! Error code


! These will eventually be provided by the host model (are not currently used,  but  will be)

real(r8), parameter :: temperature = 273
real(r8), parameter :: pressure = 100000
real(r8), parameter :: mass_density = 1

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
allocate (absTol(nSpecies))
allocate (relTol(nSpecies))
allocate (k_rateConst(nkReact))

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
! Only called at beginnning
  call k_rateConst_run(k_rateConst)

  call  chemSolve_run(nkReact, vmr_init, timeStepSize, k_rateConst, absTol, relTol, vmr_curr, ierr)

  print *, 'final_value'
  print *, vmr_curr

!-----------------------------------
! some of these will be deallocated by CPF
!-----------------------------------

  deallocate (vmr_init)
  deallocate (vmr_curr)
  deallocate (absTol)
  deallocate (relTol)
  deallocate(k_rateConst)

end program host_model_simulator
