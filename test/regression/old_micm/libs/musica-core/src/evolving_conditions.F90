! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_evolving_conditions module

!> The evolving_conditions_t type and related functions
module musica_evolving_conditions

  use musica_constants,                only : musica_dk, musica_ik
  use musica_input_output_processor,   only : input_output_processor_ptr

  implicit none
  private

  public :: evolving_conditions_t

  !> Evolving model conditions
  type evolving_conditions_t
    private
    !> Input files
    class(input_output_processor_ptr), allocatable :: input_files_(:)
  contains
    !> Get suggested output times [s]
    procedure :: get_update_times__s
    !> Update the domain state
    procedure :: update_state
    !> Preprocess the evolving conditions input data
    procedure :: preprocess_input
  end type evolving_conditions_t

  !> Constructor for evolving conditions
  interface evolving_conditions_t
    module procedure :: constructor
  end interface evolving_conditions_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for evolving conditions
  function constructor( config, domain ) result( new_obj )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_input_output_processor, only : input_output_processor_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t

    !> New evolving conditions object
    type(evolving_conditions_t), pointer :: new_obj
    !> Evolving conditions configuration data
    type(config_t), intent(inout) :: config
    !> Model domain
    class(domain_t), intent(inout) :: domain

    character(len=*), parameter :: my_name = 'evolving conditions constructor'
    type(config_t) :: file_config
    class(iterator_t), pointer :: file_iter
    type(string_t) :: file_name, file_type
    type(string_t), allocatable :: str_array(:)
    integer :: i_file, n_files

    allocate( new_obj )

    ! load the evolving conditions files
    n_files = 0
    file_iter => config%get_iterator( )
    do while( file_iter%next( ) )
      n_files = n_files + 1
    end do
    allocate( new_obj%input_files_( n_files ) )
    call file_iter%reset( )
    i_file = 1
    do while( file_iter%next( ) )
      call config%get( file_iter, file_config, my_name )
      file_name = config%key( file_iter )
      str_array = file_name%split( "." )
      file_type = str_array( size( str_array ) )%to_lower( )
      call file_config%add( "type", file_type, my_name )
      call file_config%add( "intent", "input", my_name )
      call file_config%add( "default variability", "tethered", my_name )
      call file_config%add( "file name", file_name, my_name )
      new_obj%input_files_( i_file )%val_ =>                                  &
          input_output_processor_t( file_config, domain )
      i_file = i_file + 1
    end do

    deallocate( file_iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get suggested update times [s]
  !!
  !! These times correspond to the entries in the input data
  !!
  function get_update_times__s( this ) result( times )

    use musica_array,                  only : merge_series

    !> Suggested update times
    real(kind=musica_dk), allocatable :: times(:)
    !> Evolving conditions
    class(evolving_conditions_t), intent(inout) :: this

    integer :: i_file

    allocate( times( 0 ) )
    do i_file = 1, size( this%input_files_ )
      times = merge_series( times,                                            &
                           this%input_files_( i_file )%val_%entry_times__s( ) )
    end do

  end function get_update_times__s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the model state for a given time
  subroutine update_state( this, domain, state, time__s )

    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t

    !> Evolving conditions
    class(evolving_conditions_t), intent(inout) :: this
    !> Model domain
    class(domain_t), intent(in) :: domain
    !> Domain state to update
    class(domain_state_t), intent(inout) :: state
    !> Current simulation time [s]
    real(kind=musica_dk), intent(in) :: time__s

    integer :: i_file

    do i_file = 1, size( this%input_files_ )
      call this%input_files_( i_file )%val_%update_state( domain, state,      &
                                                          time__s )
    end do

  end subroutine update_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Preprocess evolving conditions input data
  subroutine preprocess_input( this, config, domain, start_time__s,           &
      end_time__s, output_path )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t
    use musica_input_output_processor, only : input_output_processor_t
    use musica_string,                 only : string_t, to_char

    !> Evolving conditions
    class(evolving_conditions_t), intent(in) :: this
    !> Evolving conditions configuration data
    type(config_t), intent(out) :: config
    !> Model domain
    class(domain_t), intent(inout) :: domain
    !> Simulation start time [s]
    real(kind=musica_dk), intent(in) :: start_time__s
    !> Simulation end time [s]
    real(kind=musica_dk), intent(in) :: end_time__s
    !> Folder to save preprocessed data to
    character(len=*), intent(in) :: output_path

    character(len=*), parameter :: my_name = "Evolving conditions preprocessor"
    integer(kind=musica_ik) :: i_file, i_var, i_time
    type(config_t) :: config_file, file_opts
    type(string_t) :: units, file_name
    type(string_t), allocatable :: var_names(:)
    class(domain_state_t), pointer :: state
    class(input_output_processor_t), pointer :: evo_cond_file
    real(kind=musica_dk), allocatable :: times(:)

    write(*,*) "Saving evolving conditions..."

    call assert( 222456605, allocated( this%input_files_ ) )
    if( size( this%input_files_ ) .eq. 0 ) return
    do i_file = 1, size( this%input_files_ )
      var_names = this%input_files_( i_file )%val_%musica_variable_names( )
      if( size( var_names ) .eq. 0 ) cycle
      file_name = "evolving_conditions_"//trim( to_char( i_file ) )//".nc"
      write(*,*) "  - saving file '"//file_name%to_char( )//"'"
      call this%input_files_( i_file )%val_%preprocess_input( config_file,    &
                                                              "tethered" )
      call config%add( file_name%to_char( ), config_file, my_name )
      call file_opts%empty( )
      call file_opts%add( "intent", "output", my_name )
      call file_opts%add( "type", "netcdf", my_name )
      call file_opts%add( "file name", output_path//file_name%to_char( ),     &
                          my_name )
      evo_cond_file => input_output_processor_t( file_opts )
      do i_var = 1, size( var_names )
        associate( var_name => var_names( i_var ) )
        units = domain%units( var_name%to_char( ) )
        call evo_cond_file%register_output_variable( domain,                  & ! - model domain
                                                     var_name%to_char( ),     & ! - variable name
                                                     units%to_char( ) )         ! - units
        end associate
      end do
      times = this%input_files_( i_file )%val_%entry_times__s( )
      do i_time = 1, size( times )
      associate( time => times( i_time ) )
        if( time .lt. start_time__s .or. time .gt. end_time__s ) cycle
        state => domain%new_state( )
        call this%input_files_( i_file )%val_%update_state( domain, state,    &
                                                            time )
        call evo_cond_file%output( time, domain, state )
        deallocate( state )
      end associate
      end do
      deallocate( evo_cond_file )
    end do

    write(*,*) "... done!"

  end subroutine preprocess_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_evolving_conditions
