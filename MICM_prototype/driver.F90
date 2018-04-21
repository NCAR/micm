PROGRAM chemistry_cpf_mock_for_box_time_step_advance

use chemistry_specification, only: n_prognostic_variables
use solver_specification, only: nStep => ndiv
use rate_constants, only: compute_rate_constants
use ode_solve, only:  time_advance

use precision, only : r8

real(r8) :: vmr(n_prognostic_variables)
real(r8) :: advanced_vmr(n_prognostic_variables)

! this system is unstable, so 3e+1 fails
 real(r8) :: time_step_size = 2.e1_r8 ! seconds
!real(r8) :: time_step_size = .01_r8 ! seconds

! convergence criteria will have to be set somewhere and passed to ode solver.
real(r8) :: Tstart, Tend, Time
real(r8) :: RelTol(n_prognostic_variables)
real(r8) :: AbsTol(n_prognostic_variables)
integer  :: ode_return_value


! data from cpf

  vmr(:) = (/ 1._r8, 0._r8, 0._r8 /)
  print *, 'initial value'
  print *, vmr

! Only called at beginnning
    call compute_rate_constants()

    Tstart = 0._r8
    Tend   = Tstart + real(nStep,r8) * time_step_size
    Time   = Tstart
    AbsTol(:) = 1.e-9_r8
    RelTol(:) = 1.e-4_r8
! Called once to advance vmr for time_step_size seconds
    do
      if( Time < Tend ) then
        call time_advance( vmr, time_step_size, advanced_vmr, &
                           AbsTol, RelTol, ode_return_value )
        if( ode_return_value /= 1 ) then
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

END PROGRAM chemistry_cpf_mock_for_box_time_step_advance
