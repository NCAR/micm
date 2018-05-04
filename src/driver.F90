program micm_driver

!--------------------------------------------------------------------------------
! Program which prototypes the MICM driver
!--------------------------------------------------------------------------------
use external_fields,     only: set_externals
use chem_solve,          only: chem_solve_register,  chem_solve_run
use k_rate_const_module, only: k_rate_const_register, k_rate_const_init, k_rate_const_run

! This probably needs to be replaced
use solver_specification, only: nStep => ndiv


use precision, only : r8

implicit none

integer :: nSpecies                    ! Number of chemical species in run
integer :: nkReact                     ! Number of k reactions in run
real(r8), pointer :: k_rate_const(:)

real(r8),pointer :: vmr(:)
real(r8),pointer :: advanced_vmr(:)

! this system is unstable, so 3e+1 fails
real(r8) :: time_step_size = 2e+1 ! seconds

! convergence criteria will have to be set somewhere and passed to ode solver.
real(r8)         :: Tstart, Tend, Time
real(r8),pointer :: RelTol(:)
real(r8),pointer :: AbsTol(:)
integer :: ode_retcode

!-----------------------------------------------
! Initialize the chemistry packages
!-----------------------------------------------

call chem_solve_register(nSpecies)

call k_rate_const_register(nkReact)
call k_rate_const_init(nkReact, k_rate_const)


!-----------------------------------------------
! Allocate the local variables  (This will be done via CPF?)
!-----------------------------------------------
allocate (vmr(nSpecies))
allocate (advanced_vmr(nSpecies))
allocate (AbsTol(nSpecies))
allocate (RelTol(nSpecies))
allocate(k_rate_const(nkReact))


!-----------------------------------------------
! Explicitly specify the data which will come from the  cpf
!-----------------------------------------------

  vmr(:) =(/1._r8, 0._r8, 0._r8/)
  print *, 'initial value'
  print *, vmr


!-----------------------------------------------
! Simulate the XML file which CCPP will use to drive the model
!-----------------------------------------------
! Only called at beginnning
    call k_rate_const_run(k_rate_const)

  Tstart = 0._r8
  Tend   = Tstart + real(nstep,r8) * time_step_size
  Time   = Tstart
  AbsTol(:) = 1.e-9_r8
  RelTol(:) = 1.e-4_r8

! Called once to advance vmr for time_step_size seconds
    do
      if( Time < Tend ) then
        call  chem_solve_run(nkReact, vmr, time_step_size, k_rate_const, AbsTol, RelTol, advanced_vmr, ode_retcode)

        if( ode_retcode /= 1 ) then
          exit
        endif
        vmr(:) = advanced_vmr(:)
      else
        exit
      endif
      Time = Time + time_step_size
    end do


! data back to cpf

  print *, 'final_value'
  print *, advanced_vmr

!-----------------------------------
! deallocate variables
!-----------------------------------

  deallocate (vmr, advanced_vmr, k_rate_const)


end program micm_driver
