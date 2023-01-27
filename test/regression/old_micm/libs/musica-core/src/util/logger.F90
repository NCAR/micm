! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_logger module

!> The logger_t type and related functions
module musica_logger

  use musica_constants,                only : musica_dk, musica_ik

  implicit none
  private

  public :: logger_t

  !> Minimum fraction of simulation time to output progress at
  real(kind=musica_dk), parameter :: kOutputMin = 0.1_musica_dk
  !> Maximum fraction of simulation time to output progress at
  real(kind=musica_dk), parameter :: kOutputMax = 0.001_musica_dk
  !> Base output time [s]
  real(kind=musica_dk), parameter :: kBaseOutputInterval__s = 10.0_musica_dk

  !> Logger of model operational details
  type :: logger_t
    private
    !> File unit to output log to
    integer(kind=musica_ik) :: file_unit_ = 6
    !> Simulation start time [s]
    real(kind=musica_dk) :: simulation_start__s_
    !> Simulation end time [s]
    real(kind=musica_dk) :: simulation_end__s_
    !> Computation start time [s]
    real(kind=musica_dk) :: computation_start__s_
    !> Next simulation time to output at
    real(kind=musica_dk) :: next_output__s_
    !> Whether to write the progress header on the next progress output
    logical :: write_progress_header_ = .true.
  contains
    !> Log model run time and estimated time remaining
    procedure :: progress
  end type

  !> Constructor
  interface logger_t
    module procedure :: constructor
  end interface logger_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of logger_t objects
  function constructor( simulation_start__s, simulation_end__s, file_unit )   &
      result( new_obj )

    !> New logger_t object
    type(logger_t) :: new_obj
    !> Simulation start time [s]
    real(kind=musica_dk), intent(in) :: simulation_start__s
    !> Simulation end time [s]
    real(kind=musica_dk), intent(in) :: simulation_end__s
    !> File unit to output log to
    integer(kind=musica_ik), intent(in), optional :: file_unit

    new_obj%simulation_start__s_  = simulation_start__s
    new_obj%simulation_end__s_    = simulation_end__s
    new_obj%next_output__s_       = simulation_start__s
    if( present( file_unit ) ) new_obj%file_unit_ = file_unit
    call cpu_time( new_obj%computation_start__s_ )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Log the progress of the simulation and estimate time remaining
  subroutine progress( this, simulation_time__s )

    !> Logger
    class(logger_t), intent(inout) :: this
    !> Current simulation time [s]
    real(kind=musica_dk), intent(in) :: simulation_time__s

    character(len=*), parameter :: fmt = "(I11,7X,I11,6X,I11,5X, I11)"
    real(kind=musica_dk) :: sim_prog, sim_left, comp_time, comp_left, next_out
    integer(kind=musica_ik) :: f

    if( simulation_time__s .lt. this%next_output__s_ ) return
    f = this%file_unit_
    if( this%write_progress_header_ ) then
      write(f,*) "|  Simulation  |     Simulation     | Computation | Estimated time |"
      write(f,*) "| progress [s] | time remaining [s] |  time [ms]  | remaining [ms] |"
      this%write_progress_header_ = .false.
    end if
    sim_prog = simulation_time__s - this%simulation_start__s_
    sim_left = this%simulation_end__s_ - simulation_time__s
    call cpu_time( comp_time )
    comp_time = comp_time - this%computation_start__s_
    comp_left = comp_time / ( ( sim_prog + sim_left ) / sim_left - 1.0 )
    if( sim_prog .eq. 0.0 ) comp_left = 0
    write(f,fmt) int(sim_prog), int(sim_left),                                &
                 int(comp_time*1000.0_musica_dk),                             &
                 int(comp_left*1000.0_musica_dk)
    next_out = kBaseOutputInterval__s / comp_left * sim_left
    next_out = min( next_out, ( sim_prog + sim_left ) * kOutputMin )
    next_out = max( next_out, ( sim_prog + sim_left ) * kOutputMax )
    this%next_output__s_ = min( simulation_time__s + next_out,                &
                                this%simulation_end__s_ )

  end subroutine progress


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_logger

