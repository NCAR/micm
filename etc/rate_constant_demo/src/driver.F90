program rate_constant_demo

  use micm_kinetics,                   only : kinetics_t
  use micm_environment,                only : environment_t
  use musica_config,                   only : config_t

  !> Kinetics
  class(kinetics_t), pointer :: kinetics
  !> Environmental state
  type(environment_t) :: environment
  !> Path to the configuration file
  character(len=256) :: config_file_name
  !> Configuration data
  type(config_t) :: config

  ! Get the model configuration file from the command line
  if( command_argument_count( ) .ne. 1 ) then
    write(*,*) "Usage: ./rate_constant_demo configuration_file.json"
    stop 3
  end if
  call get_command_argument( 1, config_file_name )
  call config%from_file( config_file_name )

  ! do all the initialization stuff

  ! initialize the kinetics_t object
  kinetics => kinetics_t( config )

  ! do time looping

    ! get new environmental conditions from the host model
    environment%temperature = 274.5
    environment%pressure    = 100428.4

    ! call the function that calculates rate constants
    call kinetics%update_for_new_environmental_state( environment )

  ! end time loop

  ! finalize MICM

end program rate_constant_demo
