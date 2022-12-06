! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> A driver for performance evaluation of MICM
program performance_test

  use musica_constants,                only : dk => musica_dk
  use factor_solve_utilities,          only : number_of_species
  use kinetics_utilities,              only : number_of_photolysis_reactions
  use constants,                       only : kNumberOfGridCells, STREAM0, VLEN

  implicit none

  integer, parameter :: kNUmberOfTimeSteps = 100
  real(kind=dk), parameter :: kTimeStep__min = 5.0_dk

  call run_test( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Runs MICM under prescribed conditions and collect performance metrics
  subroutine run_test( )

    use micm_environment,              only : environment_t
    use micm_kinetics,                 only : kinetics_t
    use micm_ODE_solver,               only : ODE_solver_t
    use micm_ODE_solver_factory,       only : ODE_solver_builder
    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t, to_char
    use omp_lib

    character(len=*), parameter :: my_name = "MICM performance test"
    type(config_t) :: solver_config
    class(ODE_solver_t), pointer :: solver => null( )
    class(kinetics_t), pointer :: kinetics => null( )
    type(string_t), allocatable :: species_names(:)
    type(string_t), allocatable :: reaction_names(:)
    type(string_t), allocatable :: photo_reaction_names(:)
    ! Species number densities [molecule cm-3] (grid cell, species)
    real(kind=dk), allocatable :: number_densities__molec_cm3(:,:)
    type(environment_t) :: env(kNumberOfGridCells)
    integer :: i_time, error_flag, nspecies
    real(kind=dk) :: t_start, t_end

    ! Set up the kinetics calculator
    kinetics => kinetics_t( )
    call kinetics%species_names( species_names )
    call kinetics%reaction_names( reaction_names )
    call kinetics%photolysis_reaction_names( photo_reaction_names )

    ! Set up the solver for the test
    call solver_config%empty( )
    call solver_config%add( "type", "Rosenbrock", my_name )
    call solver_config%add( "chemistry time step", "min", kTimeStep__min, &
                            my_name )
    call solver_config%add( "absolute tolerance", 1.0e-12, my_name )
    call solver_config%add( "relative tolerance", 1.0e-4, my_name )
    call solver_config%add( "number of variables", number_of_species, &
                            my_name )

    ! TODO determine if the solver needs to know the number of grid cells
    !      during initialization
    solver => ODE_solver_builder( solver_config )

    ! Set up the state data strutures
    allocate( number_densities__molec_cm3( kNumberOfGridCells, &
                                           number_of_species ) )
    !$acc enter data create(number_densities__molec_cm3,env) async(STREAM0)

    ! Solve chemistry for each grid cell and time step
    do i_time = 1, kNumberOfTimeSteps
      t_start = omp_get_wtime()
      call update_species( number_densities__molec_cm3, species_names, i_time )
      call update_environment( env, photo_reaction_names, i_time )
      call kinetics%update( env )
      call solver%solve( TStart      = 0.0_dk,                                &
                         TEnd        = kTimeStep__min,                        &
                         y           = number_densities__molec_cm3,           &
                         theKinetics = kinetics,                              &
                         IErr        = error_flag )
      t_end = omp_get_wtime()
      call assert_msg( 366068772, error_flag == 0,                            &
                       "Chemistry solver failed with code "//                 &
                       to_char( error_flag ) )
      call output_state( number_densities__molec_cm3, species_names, i_time )
      write(*,*), "solve time", t_end - t_start
    end do

    !$acc exit data delete(number_densities__molec_cm3) async(STREAM0)

  end subroutine run_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the species concentrations for a given time step
  subroutine update_species( number_densities__molec_cm3, species_names,      &
      time_step )

    use musica_string,                 only : string_t

    !> Species number densities [molecule cm-3] (grid cell, species)
    real(kind=dk),  intent(inout) :: number_densities__molec_cm3(kNumberOfGridCells,number_of_species)
    type(string_t), intent(in)    :: species_names(number_of_species)
    integer,        intent(in)    :: time_step

    ! Local variables
    integer :: i, k
    real(kind=dk) :: perturb

    ! TODO determine how we want to set species concentrations

    perturb = 1._dk * time_step / ( 1._dk * kNUmberOfTimeSteps )
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector collapse(2)
    do k = 1, number_of_species 
       do i = 1, kNumberOfGridCells
          number_densities__molec_cm3(i,k) = 1.0e-6_dk * perturb
       end do
    end do
    !$acc end parallel

  end subroutine update_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the temperature, pressure, and photolysis reaction rate constants
  !! for a given time step
  subroutine update_environment( environment, photo_reaction_names, time_step )

    use micm_environment,              only : environment_t
    use musica_string,                 only : string_t

    type(environment_t), intent(inout) :: environment(kNumberOfGridCells)
    type(string_t),      intent(in)    :: photo_reaction_names(number_of_photolysis_reactions)
    integer,             intent(in)    :: time_step

    integer :: i_env, i
    real(kind=dk) :: perturb

    ! TODO determine how we want to set photolysis reaction rate constants

    perturb = 1._dk * time_step / ( 1._dk * kNUmberOfTimeSteps )
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector
    do i_env = 1, kNumberOfGridCells
       environment( i_env )%temperature = 298.15_dk * perturb ! [K]
       environment( i_env )%pressure    = 101325.0_dk * perturb ! [Pa]
       do i = 1, number_of_photolysis_reactions
          environment( i_env )%photolysis_rate_constants(i) = 1.0e-3_dk * perturb ! [s-1]
       end do
    end do
    !$acc end parallel

  end subroutine update_environment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output the state at a given time step
  subroutine output_state( number_densities__molec_cm3, species_names,        &
      time_step )

    use musica_string,                 only : string_t

    !> Species number densities [molecule cm-3] (grid cell, species)
    real(kind=dk),  intent(inout) :: number_densities__molec_cm3(:,:)
    type(string_t), intent(in)    :: species_names(:)
    integer,        intent(in)    :: time_step

    integer :: i_species

    ! TODO determine if/how we want to output state data
    write(*,*) "time step", time_step*kTimeStep__min
    do i_species = 1, size( species_names )
      if (size(number_densities__molec_cm3,1) < 100) then
          write(*,*) species_names( i_species ),                              &
                     number_densities__molec_cm3(:,i_species)
      else
          write(*,*) species_names( i_species ),                              &
                     number_densities__molec_cm3(1:100,i_species)
      end if
    end do

  end subroutine output_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program performance_test
