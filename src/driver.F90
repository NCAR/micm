program micm_driver

!--------------------------------------------------------------------------------
! Program which prototypes the MICM driver
!--------------------------------------------------------------------------------
use external_fields,     only: set_externals
use chem_solve,          only: chem_solve_register,  chem_solve_run
use k_rateConst_module, only: k_rateConst_register, k_rateConst_init, k_rateConst_run

! This probably needs to be replaced
use solver_specification, only: nStep => ndiv


use precision, only : r8

implicit none

integer :: nSpecies                    ! Number of chemical species in run
integer :: nkReact                     ! Number of k reactions in run
real(r8), pointer :: k_rateConst(:)

real(r8),pointer :: vmr(:)
real(r8),pointer :: advanced_vmr(:)

! this system is unstable, so 3e+1 fails
real(r8) :: time_step_size = 2e+1_r8 ! seconds

! convergence criteria will have to be set somewhere and passed to ode solver.
real(r8)         :: Tstart, Tend, Time
real(r8),pointer :: RelTol(:)
real(r8),pointer :: AbsTol(:)
integer :: ode_retcode

!-----------------------------------------------
! Register the chemistry packages
!-----------------------------------------------

call chem_solve_register(nSpecies)
call k_rateConst_register(nkReact)


!-----------------------------------------------
! Allocate the local variables  (This will be done via CPF?)
!-----------------------------------------------
allocate (vmr(nSpecies))
allocate (advanced_vmr(nSpecies))
allocate (AbsTol(nSpecies))
allocate (RelTol(nSpecies))
allocate(k_rateConst(nkReact))

!-----------------------------------------------
! Initialize the chemistry packages
!-----------------------------------------------

call k_rateConst_init(nkReact, k_rateConst)

!-----------------------------------------------
! Explicitly specify the data which will come from the  cpf
!-----------------------------------------------

  vmr(:) =(/1._r8, 0._r8, 0._r8/)
  print *, 'initial value'
  print *, vmr


!-----------------------------------------------
! Simulate the XML file which CCPP will use to drive the model
! Only called at beginnning
!-----------------------------------------------
    call k_rateConst_run(k_rateConst)

  Tstart = 0._r8
  Tend   = Tstart + real(nstep,r8) * time_step_size
  Time   = Tstart
  AbsTol(:) = 1.e-9_r8
  RelTol(:) = 1.e-4_r8

! Called once to advance vmr for time_step_size seconds
    do
      if( Time < Tend ) then
        call  chem_solve_run(nkReact, vmr, time_step_size, k_rateConst, AbsTol, RelTol, advanced_vmr, ode_retcode)

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
! some of these will be deallocated by CPF
!-----------------------------------

  deallocate (vmr)
  deallocate (advanced_vmr)
  deallocate (AbsTol)
  deallocate (RelTol)
  deallocate(k_rateConst)


end program micm_driver
