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
  use constants,                       only : kBoltzmann, dfactor, &
                                              ncell=>kNumberOfGridCells, &
                                              btime, etime, dtime, &
                                              blat, elat, dlat, &
                                              blon, elon, dlon, &
                                              blev, elev, dlev, &
                                              nlon, nlat, nlev, ntime, &
                                              STREAM0, VLEN, masterproc
  use netcdf
#ifdef USE_MPI
  use mpi
#endif
#else
  use constants,                       only : ncell=>kNumberOfGridCells, &
                                              STREAM0, VLEN, masterproc
#endif

  implicit none

  integer, parameter :: kNUmberOfTimeSteps = 30 
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
    type(environment_t) :: env(ncell)
    integer :: i_time, error_flag, nspecies
    real(kind=dk) :: t_start, t_end, t_measure, t_max
    integer :: myrank, mpisize, stride
    integer :: i, j, k, m, n
#ifdef USE_NETCDF
    integer :: ind, ncid
    integer :: varid(number_of_species)
    ! netCDF file name to store MICM output
    character(len=128) :: file_name = "./test_output.nc"
    real(kind=rk) :: cam_vars(ncell,dtime,number_of_species)
    real(kind=rk) :: cam_photolysis_rates(ncell,dtime,number_of_photolysis_reactions)    
    real(kind=rk), dimension(ncell,dtime) :: cam_pmid, cam_temp
#endif
#if (defined USE_NETCDF && defined USE_MPI)
    integer :: ierror

    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)
    call mpi_comm_size(mpi_comm_world, mpisize, ierror)
#else
    myrank = masterproc
    mpisize = 1
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
    ! 1. Read in the input data from CAM netCDF output
    ! 2. Replicate the data if necessary
    call read_CAM_output( cam_vars, cam_photolysis_rates, &
                          cam_temp, cam_pmid, species_names )
    if ( myrank == masterproc ) then
       ! Generate an empty netCDF file for output
       call create_netcdf_file ( file_name, species_names, ncid, varid )
    end if
#endif

    ! Set up the state data strutures
    allocate( number_densities__molec_cm3(ncell, number_of_species) )
    !$acc enter data create(number_densities__molec_cm3,env) async(STREAM0)
#ifdef USE_NETCDF
    !$acc enter data copyin(cam_vars,cam_photolysis_rates, &
    !$acc                   cam_temp,cam_pmid) async(STREAM0) 
#endif

    ! Solve chemistry for each grid cell and time step
    do i_time = 1, kNumberOfTimeSteps
      t_start = omp_get_wtime()
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
      t_end = omp_get_wtime()
      call assert_msg( 366068772, error_flag == 0,                            &
                       "Chemistry solver failed with code "//                 &
                       to_char( error_flag ) )
      if ( myrank == masterproc ) then
#ifdef USE_NETCDF
         call output_state( number_densities__molec_cm3, species_names, i_time,  &
                         ncid, varid )
#else
         call output_state( number_densities__molec_cm3, species_names, i_time )
#endif
      end if
      t_measure = t_end - t_start
#ifdef USE_MPI
      call mpi_reduce(t_measure, t_max, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, masterproc, mpi_comm_world, ierror)
#else
      t_max = t_measure
#endif
      if ( myrank == masterproc ) write(*,*) "Calculated grids per second: ", mpisize*ncell/t_max
    end do

#ifdef USE_NETCDF
    !$acc exit data delete (cam_vars,cam_photolysis_rates, &
    !$acc                   cam_temp,cam_pmid) async(STREAM0) 
#endif
    !$acc exit data delete(number_densities__molec_cm3,env) async(STREAM0)
    deallocate(number_densities__molec_cm3)

#ifdef USE_NETCDF
    ! Close the netCDF file. This frees up any internal netCDF resources
    ! associated with the file, and flushes any buffers.
    if ( myrank == masterproc ) then
       call check( nf90_close(ncid) )
       write(*,*) "Successfully write MICM output to ", file_name
    end if
#ifdef USE_MPI
    call mpi_finalize(ierror)
    if (ierror /= MPI_SUCCESS) write(*,*) "Failed to finalize MPI..."
#endif
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

    integer,          intent(in)    :: time_step
    !> Species number densities [molecule cm-3] (grid cell, species)
    real(kind=dk),    intent(inout) :: number_densities__molec_cm3(ncell,number_of_species)
#ifdef USE_NETCDF
    !> Spevies volume mixing ratios from CAM [mol mol-1]
    real(kind=rk),       intent(in) :: cam_vars(ncell,dtime,number_of_species)
    type(environment_t), intent(in) :: environment(ncell)
#endif

    !> Local variables
#ifdef USE_NETCDF
    integer :: n, m, p
#else
    integer :: i, k
    real(kind=dk) :: perturb
#endif

    ! TODO determine how we want to set species concentrations

#ifdef USE_NETCDF
    p = mod(time_step, dtime)
    if (p == 0) p = p + dtime
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector collapse(2)
    do m = 1, number_of_species
       do n = 1, ncell
          !> Convert unit from mol/mol to molecule/cm3
          number_densities__molec_cm3(n,m) &
                = real(cam_vars(n,p,m), kind=dk) * 10._dk * &
                  environment(n)%pressure / &
                  (kBoltzmann*1.e7_dk*environment(n)%temperature)
       end do
    end do
    !$acc end parallel
#else
    perturb = 1._dk * time_step / ( 1._dk * kNUmberOfTimeSteps )
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector collapse(2)
    do k = 1, number_of_species 
       do i = 1, ncell 
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
  subroutine update_environment( environment, time_step, cam_photolysis_rates, &
                                 cam_temp, cam_pmid )
#else
  subroutine update_environment( environment, time_step )
#endif

    use micm_environment,              only : environment_t

    integer,             intent(in)    :: time_step
    type(environment_t), intent(inout) :: environment(ncell)
#ifdef USE_NETCDF
    real(kind=rk),       intent(in)    :: cam_photolysis_rates(ncell,dtime,number_of_photolysis_reactions)
    real(kind=rk),       intent(in)    :: cam_temp(ncell,dtime)
    real(kind=rk),       intent(in)    :: cam_pmid(ncell,dtime)
#endif

    !> Local variables
#ifdef USE_NETCDF
    integer :: i_env, m, p
#else
    integer :: i_env, i
    real(kind=dk) :: perturb
#endif
    ! TODO determine how we want to set photolysis reaction rate constants

#ifdef USE_NETCDF
    p = mod(time_step, dtime)
    if (p == 0) p = p + dtime
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector
    do i_env = 1, ncell
       environment( i_env )%temperature = real(cam_temp(i_env,p), kind=dk) ! [K]
       environment( i_env )%pressure    = real(cam_pmid(i_env,p), kind=dk) ! [Pa]
       do m = 1, number_of_photolysis_reactions
          environment( i_env )%photolysis_rate_constants(m) &
                          = real(cam_photolysis_rates(i_env,p,m), kind=dk) ! [s-1]
       end do
    end do
    !$acc end parallel
#else
    perturb = 1._dk * time_step / ( 1._dk * kNUmberOfTimeSteps )
    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector
    do i_env = 1, ncell
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
    real(kind=dk),  intent(in) :: number_densities__molec_cm3(ncell,number_of_species)
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
      if (ncell < 10) then
          write(*,*) species_names( i_species ),                              &
                     number_densities__molec_cm3(:,i_species)
      else
          write(*,*) species_names( i_species ),                              &
                     number_densities__molec_cm3(1:10,i_species)
      end if
    end do

#ifdef USE_NETCDF
    !> These settings tell netcdf to write which timestep, as determined by start(2)
    counts = (/int(ncell/dfactor), 1/)
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
  real(kind=rk), intent(out) :: cam_vars(ncell,dtime,number_of_species)
  !> Photolysis reaction rates from CAM (unit: 1/s)
  real(kind=rk), intent(out) :: cam_photolysis_rates(ncell,dtime,number_of_photolysis_reactions)
  !> Temperature from CAM (unit: K)
  real(kind=rk), intent(out) :: cam_temp(ncell,dtime)
  !> Pressure at mid-point of layers from CAM (unit: Pa)
  real(kind=rk), intent(out) :: cam_pmid(ncell,dtime)

  !> Local variables
  integer            :: ncid, varid, nvar
  integer            :: i, m
  character(len=256) :: file_name
  character(len=128) :: str
  logical            :: missing_var_flag
  !> Corresponding MICM photolysis reaction names in CAM output for Chapman mechanism
!  character(len=128), parameter :: cam_photo_reaction_names(number_of_photolysis_reactions) = &
!                                   (/'jo2_a','jo3_a','jo3_b'/)

!  !> Corresponding MICM photolysis reaction names in CAM output for TS1 mechanism
  character(len=*), parameter :: cam_photo_reaction_names(number_of_photolysis_reactions) = &
             (/ 'jterpnit       ', 'jo2_b          ', 'jmvk           ', 'jmgly          ', &
                'jxylenooh      ', 'jso3           ', 'jchbr3         ', 'jterprd2       ', &
                'jch4_b         ', 'jmek           ', 'jacet          ', 'jterp2ooh      ', &
                'jbrono2_a      ', 'jbigald4       ', 'jcf2clbr       ', 'jh2o_c         ', &
                'jbro           ', 'jtepomuc       ', 'jxylolooh      ', 'jso2           ', &
                'jmekooh        ', 'jch2o_b        ', 'jbrono2_b      ', 'jno            ', &
                'jphenooh       ', 'jh2o_b         ', 'jc3h7ooh       ', 'jccl4          ', &
                'jnoa           ', 'jno2           ', 'jbigald3       ', 'jch2br2        ', &
                'jbepomuc       ', 'jmacr_b        ', 'jnterpooh      ', 'jn2o           ', &
                'jrooh          ', 'jch2o_a        ', 'jhbr           ', 'jpooh          ', &
                'jhcfc22        ', 'jnc4cho        ', 'jhonitr        ', 'jxooh          ', &
                'jch4_a         ', 'jhobr          ', 'jno3_a         ', 'jh2o_a         ', &
                'jch3co3h       ', 'jhcfc141b      ', 'jhyac          ', 'jonitr         ', &
                'jterpooh       ', 'jcl2           ', 'jhcfc142b      ', 'jisopooh       ', &
                'jh2o2          ', 'jch3cho        ', 'jbigald1       ', 'joclo          ', &
                'jhno3          ', 'jc2h5ooh       ', 'jmacr_a        ', 'jco2           ', &
                'jbenzooh       ', 'jn2o5_b        ', 'jcfc115        ', 'jcfcl3         ', &
                'jterprd1       ', 'jso            ', 'jmpan          ', 'jclono2_a      ', &
                'jisopnooh      ', 'jalkooh        ', 'jcf3br         ', 'jch3cl         ', &
                'jhcl           ', 'jtolooh        ', 'jch3ccl3       ', 'jglyald        ', &
                'jeooh          ', 'jbrcl          ', 'jocs           ', 'jc6h5ooh       ', &
                'jcfc114        ', 'jho2no2_b      ', 'jpan           ', 'jhocl          ', &
                'jo3_b          ', 'jcf2cl2        ', 'jno3_b         ', 'jbigald2       ', &
                'jch3ooh        ', 'jclo           ', 'jcfc113        ', 'jho2no2_a      ', &
                'jch3br         ', 'jclono2_b      ', 'jglyoxal       ', 'jo3_a          ', &
                'jn2o5_a        ', 'jhpald         ', 'jbzooh         ', 'jh2402         ', &
                'jcl2o2         ', 'jalknit        ', 'jbigald        ' /)

  !> Path to CAM output file  
  file_name = '/glade/scratch/fvitt/archive/TS1_chem_output_t01/atm/hist/TS1_chem_output_t01.cam.h1.2010-01-06-00000.nc'

  !> Sanity check before reading in the CAM output
  call sanity_check_constants()

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
        call init_array ( cam_vars(1:ncell,1:dtime,m), 1.e-6_rk )
     else
        call get_cam_var ( ncid, varid, str, cam_vars(1:ncell,1:dtime,m) )
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
        call init_array ( cam_photolysis_rates(1:ncell,1:dtime,m), 1.e-6_rk )
     else
        call get_cam_var ( ncid, varid, str, cam_photolysis_rates(1:ncell,1:dtime,m) )
     end if
  end do

  !> Read in temperature 
  str = 'T'
  !> Get the ID of a variable, based on its name.
  call check( nf90_inq_varid(ncid, str, varid), str, missing_var_flag )
  if ( missing_var_flag ) then
     !> Miss temperature output from CAM; set to 298.15 by default
     call init_array ( cam_temp, 298.15_rk )
  else
     call get_cam_var ( ncid, varid, str, cam_temp )
  end if

  !> Read in pressure
  str = 'PMID'
  !> Get the ID of a variable, based on its name.
  call check( nf90_inq_varid(ncid, str, varid), str, missing_var_flag )
  if ( missing_var_flag ) then
     !> Miss pressure output from CAM; set to 101325 by default
     call init_array ( cam_pmid, 101325._rk )
  else
     call get_cam_var ( ncid, varid, str, cam_pmid )
  end if

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )
  write(*,*) "Successfully read in the CAM output from ", file_name

  end subroutine read_CAM_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize a multi-dimension array with a given value
  subroutine init_array ( array, val )

  real(kind=rk), dimension(ncell,dtime), intent(out) :: array
  real(kind=rk), intent(in) :: val

  ! Local variables
  integer :: i, n

  do n = 1, dtime
     do i = 1, ncell
        array(i,n) = val
     end do
  end do

  end subroutine init_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check if the constants to read in the slices of CAM output are valid
  subroutine sanity_check_constants ()

  if ( btime < 1 .or. btime > etime ) then
     write(*,*) "Invalid value for btime: ", btime
     stop 
  end if
  if ( etime > ntime ) then
     write(*,*) "Invalid value for etime: ", etime
     stop 
  end if
  if ( blev < 1 .or. blev > elev ) then
     write(*,*) "Invalid value for blev: ", blev
     stop               
  end if                
  if ( elev > nlev ) then
     write(*,*) "Invalid value for elev: ", elev
     stop
  end if
  if ( blat < 1 .or. blat > elat ) then
     write(*,*) "Invalid value for blat: ", blat
     stop
  end if
  if ( elat > nlat ) then
     write(*,*) "Invalid value for elat: ", elat
     stop
  end if
  if ( blon < 1 .or. blon > elon ) then
     write(*,*) "Invalid value for blon: ", blon
     stop
  end if
  if ( elon > nlon ) then
     write(*,*) "Invalid value for elon: ", elon
     stop
  end if

  end subroutine sanity_check_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read in the variable from CAM output and reshape to the shape for MICM
  subroutine get_cam_var ( ncid, varid, str, array )

  integer, intent(in) :: ncid, varid
  character(len=128), intent(in) :: str
  real(kind=rk), intent(out) :: array(ncell,dtime)
  
  ! Local variables
  integer, parameter :: start(4) = (/blon,blat,blev,btime/), &
                        counts(4) = (/dlon,dlat,dlev,dtime/), &
                        stride(4) = (/1,1,1,1/)
  real(kind=rk) :: var_io(dlon,dlat,dlev,dtime)
  integer :: i, n

  !> Read in the values of data
  call check( nf90_get_var(ncid, varid, var_io, start, counts, stride), str )
  !> Reshape the array
  n = ncell/dfactor
  array(1:n,1:dtime) = reshape( var_io, (/n, dtime/) )
  !> Replicate the data by "dfactor" times
  do i = 2, dfactor
     array((i-1)*n+1:i*n,1:dtime) = array(1:n,1:dtime)
  end do

  end subroutine get_cam_var

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
  call check( nf90_def_dim(ncid, "number_of_gridcells", ncell/dfactor, ncell_dimid) )

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
