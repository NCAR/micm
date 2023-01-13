! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> A driver for performance evaluation of MICM
program performance_test

  use musica_constants,                only : dk => musica_dk, &
                                              rk => musica_rk
  use factor_solve_utilities,          only : number_of_species
  use kinetics_utilities,              only : number_of_photolysis_reactions
#ifdef USE_NETCDF
  use constants,                       only : kBoltzmann, kNumberOfGridCells, &
                                              nlon, nlat, nlev, ntime, &
                                              STREAM0, VLEN
  use netcdf
#else
  use constants,                       only : kNumberOfGridCells, &
                                              STREAM0, VLEN
#endif

  implicit none

  integer, parameter :: kNUmberOfTimeSteps = 5 
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
#ifdef USE_NETCDF
    integer :: ncid
    integer :: varid(number_of_species)
    ! netCDF file name to store MICM output
    character(len=128) :: file_name = "./test_output.nc"
    real(kind=rk), dimension(:,:,:,:,:), allocatable :: cam_vars, &
                                               cam_photolysis_rates
    real(kind=rk), dimension(:,:,:,:), allocatable :: cam_pmid, cam_temp
#endif

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

#ifdef USE_NETCDF
    allocate(cam_vars(nlon,nlat,nlev,ntime,number_of_species))
    allocate(cam_photolysis_rates(nlon,nlat,nlev,ntime,number_of_photolysis_reactions))
    allocate(cam_temp(nlon,nlat,nlev,ntime))
    allocate(cam_pmid(nlon,nlat,nlev,ntime))
    ! Read in the input data from CAM netCDF output
    call read_CAM_output( cam_vars, cam_photolysis_rates, &
                          cam_temp, cam_pmid, species_names )
    ! Generate an empty netCDF file for output
    call create_netcdf_file ( file_name, species_names, ncid, varid )
#endif

    ! Set up the state data strutures
    allocate( number_densities__molec_cm3( kNumberOfGridCells, &
                                           number_of_species ) )
    !$acc enter data create(number_densities__molec_cm3,env) async(STREAM0)

    ! Solve chemistry for each grid cell and time step
    do i_time = 1, kNumberOfTimeSteps
!      t_start = omp_get_wtime()
#ifdef USE_NETCDF
      call update_environment( env, i_time, cam_photolysis_rates, cam_temp, cam_pmid )
      call update_species( number_densities__molec_cm3, i_time, cam_vars, env )
#else
      call update_environment( env, i_time )
      call update_species( number_densities__molec_cm3, i_time )
#endif
      call kinetics%update( env )
      call solver%solve( TStart      = 0.0_dk,                                &
                         TEnd        = kTimeStep__min,                        &
                         y           = number_densities__molec_cm3,           &
                         theKinetics = kinetics,                              &
                         IErr        = error_flag )
!      t_end = omp_get_wtime()
      call assert_msg( 366068772, error_flag == 0,                            &
                       "Chemistry solver failed with code "//                 &
                       to_char( error_flag ) )
#ifdef USE_NETCDF
      call output_state( number_densities__molec_cm3, species_names, i_time, ncid, varid )
#else
      call output_state( number_densities__molec_cm3, species_names, i_time )
#endif
      write(*,*) "solve time", t_end - t_start
    end do

#ifdef USE_NETCDF
    !$acc exit data delete(number_densities__molec_cm3,cam_vars, &
    !$acc                  cam_photolysis_rates,cam_temp,cam_pmid) async(STREAM0)

    ! Close the netCDF file. This frees up any internal netCDF resources
    ! associated with the file, and flushes any buffers.
    call check( nf90_close(ncid) )
    write(*,*) "Successfully write MICM output to ", file_name
#else
    !$acc exit data delete(number_densities__molec_cm3) async(STREAM0)
#endif

  end subroutine run_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the species concentrations for a given time step
#ifdef USE_NETCDF
  subroutine update_species( number_densities__molec_cm3, time_step, cam_vars, environment )
#else
  subroutine update_species( number_densities__molec_cm3, time_step )
#endif

    use micm_environment,              only : environment_t

    !> Species number densities [molecule cm-3] (grid cell, species)
    real(kind=dk),    intent(inout) :: number_densities__molec_cm3(kNumberOfGridCells,number_of_species)
    integer,          intent(in)    :: time_step
#ifdef USE_NETCDF
    !> Spevies volume mixing ratios from CAM [mol mol-1]
    real(kind=rk),       intent(in) :: cam_vars(nlon,nlat,nlev,ntime,number_of_species)
    type(environment_t), intent(in) :: environment(kNumberOfGridCells)
#endif

    !> Local variables
#ifdef USE_NETCDF
    integer :: i, j, k, m, n, p
#else
    integer :: i, k
    real(kind=dk) :: perturb
#endif

    ! TODO determine how we want to set species concentrations

#ifdef USE_NETCDF
    p = mod(time_step, ntime)
    if (p == 0) p = p + ntime
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector collapse(4)
    do m = 1, number_of_species
       do k = 1, nlev
          do j = 1, nlat
             do i = 1, nlon
                n = i + (j-1) * nlon + (k-1) * nlat * nlon
                !> Convert unit from mol/mol to molecule/cm3
                number_densities__molec_cm3(n,m) &
                      = real(cam_vars(i,j,k,p,m), kind=dk) * 10._dk * &
                        environment(n)%pressure / &
                        (kBoltzmann*1.e7_dk*environment(n)%temperature)
             end do
          end do
       end do
    end do
    !$acc end parallel
#else
    perturb = 1._dk * time_step / ( 1._dk * kNUmberOfTimeSteps )
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector collapse(2)
    do k = 1, number_of_species 
       do i = 1, kNumberOfGridCells
          number_densities__molec_cm3(i,k) = 1.0e-6_dk * perturb
       end do
    end do
    !$acc end parallel
#endif

  end subroutine update_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the temperature, pressure, and photolysis reaction rate constants
  !! for a given time step
#ifdef USE_NETCDF
  subroutine update_environment( environment, time_step, &
                                 cam_photolysis_rates, cam_temp, cam_pmid )
#else
  subroutine update_environment( environment, time_step )
#endif

    use micm_environment,              only : environment_t

    type(environment_t), intent(inout) :: environment(kNumberOfGridCells)
    integer,             intent(in)    :: time_step
#ifdef USE_NETCDF
    real(kind=rk),       intent(in)    :: cam_photolysis_rates(nlon,nlat,nlev,ntime,number_of_photolysis_reactions)
    real(kind=rk),       intent(in)    :: cam_temp(nlon,nlat,nlev,ntime)
    real(kind=rk),       intent(in)    :: cam_pmid(nlon,nlat,nlev,ntime)
#endif

    !> Local variables
#ifdef USE_NETCDF
    integer :: i_env, i, j, k, m, p
#else
    integer :: i_env, i
    real(kind=dk) :: perturb
#endif
    ! TODO determine how we want to set photolysis reaction rate constants

#ifdef USE_NETCDF
    p = mod(time_step, ntime)
    if (p == 0) p = p + ntime
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector collapse(3)
    do k = 1, nlev
       do j = 1, nlat
          do i = 1, nlon
             i_env = i + (j-1) * nlon + (k-1) * nlat * nlon
             environment( i_env )%temperature = real(cam_temp(i,j,k,p), kind=dk) ! [K]
             environment( i_env )%pressure    = real(cam_pmid(i,j,k,p), kind=dk) ! [Pa]
             do m = 1, number_of_photolysis_reactions
                environment( i_env )%photolysis_rate_constants(m) &
                                = real(cam_photolysis_rates(i,j,k,p,m), kind=dk) ! [s-1]
             end do
          end do
       end do
    end do
    !$acc end parallel
#else
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
#endif

  end subroutine update_environment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output the state at a given time step
#ifdef USE_NETCDF
  subroutine output_state( number_densities__molec_cm3, species_names, &
                           time_step, ncid, varid )
#else
  subroutine output_state( number_densities__molec_cm3, species_names, &
                           time_step )
#endif

    use musica_string,                 only : string_t

    !> Species number densities [molecule cm-3] (grid cell, species)
    real(kind=dk),  intent(in) :: number_densities__molec_cm3(kNumberOfGridCells,number_of_species)
    type(string_t), intent(in) :: species_names(number_of_species)
    integer,        intent(in) :: time_step
#ifdef USE_NETCDF
    integer,        intent(in) :: ncid
    integer,        intent(in) :: varid(number_of_species)
#endif

    !> Local variables
#ifdef USE_NETCDF
    integer, parameter :: ndims = 2
    integer :: counts(ndims), start(ndims)
#endif
    integer :: i_species

    ! TODO determine if/how we want to output state data
    write(*,*) "time step", time_step*kTimeStep__min
    do i_species = 1, number_of_species
      if (kNumberOfGridCells < 10) then
          write(*,*) species_names( i_species ),                              &
                     number_densities__molec_cm3(:,i_species)
      else
          write(*,*) species_names( i_species ),                              &
                     number_densities__molec_cm3(1:10,i_species)
      end if
    end do

#ifdef USE_NETCDF
    !> These settings tell netcdf to write which timestep, as determined by start(2)
    counts = (/kNumberOfGridCells, 1/)
    start = (/1, time_step/)

    ! Write the data to the file.
    do i_species = 1, number_of_species
       call check( nf90_put_var(ncid, varid(i_species), &
                                number_densities__molec_cm3(:,i_species), &
                                start=start, count=counts) )
    end do
#endif

  end subroutine output_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_NETCDF
  !> Read the CAM netCDF output as MICM input
  !> Example could be found at https://www.unidata.ucar.edu/software/netcdf/examples/programs/
  subroutine read_CAM_output ( cam_vars, cam_photolysis_rates, &
                               cam_temp, cam_pmid, species_names )

  use musica_string, only : string_t
 
  !> Input variables
  type(string_t), intent(in) :: species_names(number_of_species)
  !> Species volume mixing ratio from CAM (unit: mol/mol)
  real(kind=rk), intent(out) :: cam_vars(nlon,nlat,nlev,ntime,number_of_species)
  !> Photolysis reaction rates from CAM (unit: 1/s)
  real(kind=rk), intent(out) :: cam_photolysis_rates(nlon,nlat,nlev,ntime,number_of_photolysis_reactions)
  !> Temperature from CAM (unit: K)
  real(kind=rk), intent(out) :: cam_temp(nlon,nlat,nlev,ntime)
  !> Pressure at mid-point of layers from CAM (unit: Pa)
  real(kind=rk), intent(out) :: cam_pmid(nlon,nlat,nlev,ntime)

  !> Local variables
  integer            :: ncid, varid, nvar
  integer            :: i, j, k, n, m
  character(len=256) :: file_name
  character(len=128) :: str
  logical            :: missing_var_flag
  !> Corresponding MICM photolysis reaction names in CAM output for Chapman mechanism
  character(len=128), parameter :: cam_photo_reaction_names(number_of_photolysis_reactions) = &
                                   (/'jo2_a','jo3_a','jo3_b'/)

  !> Path to CAM output file  
  file_name = '/glade/scratch/fvitt/archive/TS1_chem_output_t01/atm/hist/TS1_chem_output_t01.cam.h1.2010-01-06-00000.nc'

  !> Open an existing netcdf file
  call check( nf90_open(file_name, nf90_nowrite, ncid) )

  !> Read in species concentration
  do m = 1, number_of_species
     !> Get the ID of a variable, based on its name.
     str = species_names(m)   ! this step is necessary to convert "string_t" type
                              ! to "character" type for netCDF subroutine call
     call check( nf90_inq_varid(ncid, str, varid), str, missing_var_flag )
     if ( missing_var_flag ) then
        !> Miss this variable from CAM output; set to 1.e-6 by default
        do n = 1, ntime
           do k = 1, nlev
              do j = 1, nlat
                 do i = 1, nlon
                    cam_vars(i,j,k,n,m) = 1.e-6_rk
                 end do
              end do
           end do
        end do
     else
        !> Read in the values of data
        call check( nf90_get_var(ncid, varid, cam_vars(:,:,:,:,m)), str )
     end if
  end do

  !> Read in photolysis reaction rates
  do m = 1, number_of_photolysis_reactions
     !> Get the ID of a variable, based on its name.
     str = cam_photo_reaction_names(m)   ! this step is necessary to convert "string_t" type
                                         ! to "character" type for netCDF subroutine call
     call check( nf90_inq_varid(ncid, str, varid), str, missing_var_flag )
     if ( missing_var_flag ) then
        !> Miss this photolysis reaction from CAM output; set to 1.e-6 by default
        do n = 1, ntime
           do k = 1, nlev
              do j = 1, nlat
                 do i = 1, nlon
                    cam_photolysis_rates(i,j,k,n,m) = 1.e-6_rk
                 end do
              end do
           end do
        end do
     else
        !> Read in the values of data
        call check( nf90_get_var(ncid, varid, cam_photolysis_rates(:,:,:,:,m)), str )
     end if
  end do

  !> Read in temperature 
  str = 'T'
  !> Get the ID of a variable, based on its name.
  call check( nf90_inq_varid(ncid, str, varid), str, missing_var_flag )
  if ( missing_var_flag ) then
     !> Miss temperature output from CAM; set to 298.15 by default
     do n = 1, ntime
        do k = 1, nlev
           do j = 1, nlat
              do i = 1, nlon
                 cam_temp(i,j,k,n) = 298.15_rk
              end do
           end do
        end do
     end do
  else
     !> Read in the values of data
     call check( nf90_get_var(ncid, varid, cam_temp), str )
  end if

  !> Read in pressure
  str = 'PMID'
  !> Get the ID of a variable, based on its name.
  call check( nf90_inq_varid(ncid, str, varid), str, missing_var_flag )
  if ( missing_var_flag ) then
     !> Miss pressure output from CAM; set to 101325 by default
     do n = 1, ntime
        do k = 1, nlev
           do j = 1, nlat
              do i = 1, nlon
                 cam_pmid(i,j,k,n) = 101325._rk
              end do
           end do
        end do
     end do
  else
     !> Read in the values of data
     call check( nf90_get_var(ncid, varid, cam_pmid), str )
  end if

  !$acc enter data copyin(cam_vars,cam_photolysis_rates,cam_temp,cam_pmid) async(STREAM0)

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )
  write(*,*) "Successfully read in the CAM output from ", file_name

  end subroutine read_CAM_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generate a netCDF file to store MICM output
  subroutine create_netcdf_file ( file_name, species_names, ncid, varid )

  use musica_string,                 only : string_t
  
  character(len=128), intent(in) :: file_name
  type(string_t),     intent(in) :: species_names(number_of_species)
  integer,           intent(out) :: ncid
  integer,           intent(out) :: varid(number_of_species)

  !> Local variables
  character(len=128) :: str
  character(len=128), parameter :: units = "molecule cm-3"
  integer, parameter :: ndims = 2
  integer :: ncell_dimid, time_dimid, i_species
  integer :: dimids(ndims)

  !> Create the netCDF file. The nf90_clobber parameter tells netCDF to
  !> overwrite this file, if it already exists.
  call check( nf90_create(file_name, nf90_clobber, ncid) )

  !> Define the dimensions. NetCDF will hand back an ID for each. 
  call check( nf90_def_dim(ncid, "time", nf90_unlimited, time_dimid) )
  call check( nf90_def_dim(ncid, "number_of_gridcells", kNumberOfGridCells, ncell_dimid) )

  !> The dimids array is used to pass the IDs of the dimensions of
  !> the variables. Note that in fortran arrays are stored in
  !> column-major format.
  dimids = (/ ncell_dimid, time_dimid /)

  !> Define the variables.
  do i_species = 1, number_of_species
     str = species_names(i_species)
     call check( nf90_def_var(ncid, str, nf90_double, dimids, varid(i_species)) )
     !> Assign units attributes to the netCDF variables.
     call check( nf90_put_att(ncid, varid(i_species), "units", units) )
  end do

  !> End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )

  end subroutine create_netcdf_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether a netCDF subroutine call succeeds or not
  subroutine check( istatus, var, missing_var_flag )
 
  integer, intent(in) :: istatus
  character(len=128), intent(in), optional :: var
  logical, intent(out), optional :: missing_var_flag

  if (istatus /= nf90_noerr) then
     if (present(var)) then
        print *, trim(nf90_strerror(istatus))//" : "//trim(adjustl(var))
        if (present(missing_var_flag)) then
           missing_var_flag = .true.
        end if
     else
        print *, trim(nf90_strerror(istatus))
     end if
!     stop "Stopped"
  else
     if (present(missing_var_flag)) then
        missing_var_flag = .false.
     end if
  end if
 
  end subroutine check
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program performance_test
