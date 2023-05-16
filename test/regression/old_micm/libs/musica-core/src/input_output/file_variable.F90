! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_variable module

!> The file_variable_t type and related functions
module musica_file_variable

  use musica_constants,                only : musica_dk, musica_ik
  use musica_convert,                  only : convert_t
  use musica_string,                   only : string_t

  implicit none
  private

  public :: file_variable_t, file_variable_ptr, find_variable_by_name,    &
            find_variable_by_musica_name

  !> @name Variability types
  !! @{
  integer(kind=musica_ik), parameter :: VARIABILITY_PROGNOSED = 1
  integer(kind=musica_ik), parameter :: VARIABILITY_TETHERED  = 2
  integer(kind=musica_ik), parameter :: VARIABILITY_FIXED     = 3
  !> @}

  !> A File variable
  !!
  !! The file_variable_t handles all conversions, offsetting, scaling,
  !! etc. and can be used to return sub-sets of the file data during the
  !! simulation in MUSICA units after applying any specified conversions.
  !!
  type, abstract :: file_variable_t
    private
    !> Name in the file
    type(string_t) :: name_
    !> Expected MUSICA name
    type(string_t) :: musica_name_
    !> Dimension names
    type(string_t), allocatable :: dimension_names_(:)
    !> Units for variable in file data
    type(string_t) :: units_
    !> MUSICA units for the variable
    type(string_t) :: musica_units_
    !> Converter to MUSICA units
    type(convert_t) :: converter_
    !> Type of variability specified in the configuration
    integer(kind=musica_ik) :: variability_ = VARIABILITY_PROGNOSED
    !> Scaling factor
    real(kind=musica_dk) :: scale_factor_ = 1.0_musica_dk
    !> Offset (applied to file data after scaling and before unit conversion)
    real(kind=musica_dk) :: offset_ = 0.0_musica_dk
    !> Shift (applied after unit conversion)
    real(kind=musica_dk) :: shift_ = 0.0_musica_dk
  contains
    !> Returns the name of the variable
    procedure :: name => variable_name
    !> Returns the MUSICA name for the variable
    procedure :: musica_name
    !> Returns the name of one of the variable dimensions
    procedure :: dimension_name
    !> Returns a flag indicating whether the variable is specified as being
    !! "prognosed" during the simulation
    procedure :: is_prognosed
    !> Returns a flag indicating whether the variable is specified as being
    !! "tethered" to input conditions
    procedure :: is_tethered
    !> Returns a flag indicating whether the variable is specified as being
    !! "fixed" during the simulation
    procedure :: is_fixed
    !> Returns the units used in the file for the variable
    procedure :: units
    !> Sets the MUSICA units for the variable
    procedure :: set_musica_units
    !> Scale, offset, convert units, and shift a set of data to musica values
    !! (Should only be called by extending types)
    procedure :: convert_to_musica_values
    !> Scale, offset, convert units, and shift a set of data to file values
    !! (Should only be called by extending types)
    procedure :: convert_to_file_values
    !> Gets the variable dimensions
    procedure(get_dimensions), deferred :: get_dimensions
    !> Gets a sub-set of the variable data for a specified index range
    !!
    !! Data are returned after applying conversions set up during
    !! initialization.
    !!
    procedure(get_data), deferred :: get_data
    !> Outputs data to the file for a given time step
    procedure(output), deferred :: output
    !> Returns a flag indicating whether two file_variable_t objects refer to
    !! the same file variable
    procedure(is_same_as), deferred :: is_same_as
    !> Prints the properties of the variable
    procedure :: print => do_print
    !> @name Equality comparison
    !! @{
    procedure, private :: equals
    generic :: operator(==) => equals
    procedure, private :: not_equals
    generic :: operator(/=) => not_equals
    !> @}
    !> Loads matching and conversion options specified in the configuration data
    procedure, private :: load_configured_options
    !> Does standard file -> MUSICA name conversions
    procedure, private :: do_standard_name_conversions
    !> Sets a shift in input/output data
    procedure, private :: set_shift
    !> Private constructor
    !! (Should only be called by constructors of extending types)
    procedure :: private_constructor
  end type file_variable_t

  !> Unique pointer to file_variable_t objects
  type :: file_variable_ptr
    class(file_variable_t), pointer :: val_ => null( )
  contains
    !> Finalizes the pointer
    final :: file_variable_ptr_finalize
  end type file_variable_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the variable dimensions
  function get_dimensions( this ) result( dimensions )
    use musica_file_dimension_range,   only : file_dimension_range_t
    import file_variable_t
    !> File dimensions
    type(file_dimension_range_t), allocatable :: dimensions(:)
    !> File variable
    class(file_variable_t), intent(in) :: this
  end function get_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets a sub-set of the data from the file
  !!
  !! Conversions are applied in the following order:
  !! - scaling
  !! - offsetting
  !! - conversion to MUSICA units
  !! - shifting
  !!
  subroutine get_data( this, file, range, values )
    use musica_constants,              only : musica_ik, musica_dk
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t
    import file_variable_t
    !> File variable
    class(file_variable_t), intent(in) :: this
    !> Input/Output file
    class(file_t), intent(inout) :: file
    !> Range of data to return
    !!
    !! Only one range in the array is permitted to have size > 1, so that
    !! results are always returned as a rank 1 array
    !!
    type(file_dimension_range_t), intent(in) :: range(:)
    !> Values to return
    real(kind=musica_dk), target, intent(out) :: values(:)
  end subroutine get_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs data to the file for a given timestep
  !!
  !! The state value will be converted according to the variable
  !! configuration prior to outputting it to the file.
  !!
  subroutine output( this, file, unlimited_dimension_value, indices,          &
      state_value )
    use musica_constants,              only : musica_dk
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t
    import file_variable_t
    !> File variable
    class(file_variable_t), intent(inout) :: this
    !> Output file
    class(file_t), intent(inout) :: file
    !> Unlimited dimension value
    real(kind=musica_dk), intent(in) :: unlimited_dimension_value
    !> Indices of data point to output in all non-unlimited dimensions
    class(file_dimension_range_t), intent(in) :: indices(:)
    !> Value to output
    real(kind=musica_dk), intent(in) :: state_value
  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether two file_variable_t objects refer to
  !! the same file variable
  logical elemental function is_same_as( a, b )
    import file_variable_t
    !> File variable a
    class(file_variable_t), intent(in) :: a
    !> File variable b
    class(file_variable_t), intent(in) :: b
  end function is_same_as

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the variable
  type(string_t) function variable_name( this )

    !> File variable
    class(file_variable_t), intent(in) :: this

    variable_name = this%name_

  end function variable_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the expected MUSICA name for the variable
  type(string_t) function musica_name( this )

    !> File variable
    class(file_variable_t), intent(in) :: this

    musica_name = this%musica_name_

  end function musica_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns one of the variable dimension names by dimension index
  type(string_t) function dimension_name( this, dimension_index )

    use musica_assert,                 only : assert

    !> File variable
    class(file_variable_t), intent(in) :: this
    !> Dimension index
    integer(kind=musica_ik), intent(in) :: dimension_index

    call assert( 477973132,                                                   &
                 dimension_index .le. size( this%dimension_names_ ) .and.     &
                 dimension_index .ge. 1 )
    dimension_name = this%dimension_names_( dimension_index )

  end function dimension_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether the variable is specified as being
  !! "prognosed" during the simulation
  logical function is_prognosed( this )

    !> File variable
    class(file_variable_t), intent(in) :: this

    is_prognosed = this%variability_ .eq. VARIABILITY_PROGNOSED

  end function is_prognosed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether the variable is specified as being
  !! "tethered" to input conditions
  logical function is_tethered( this )

    !> File variable
    class(file_variable_t), intent(in) :: this

    is_tethered = this%variability_ .eq. VARIABILITY_TETHERED

  end function is_tethered

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether the variable is specified as being
  !! "fixed" during the simulation
  logical function is_fixed( this )

    !> File variable
    class(file_variable_t), intent(in) :: this

    is_fixed = this%variability_ .eq. VARIABILITY_FIXED

  end function is_fixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units used in the file for the variable
  type(string_t) function units( this )

    !> File variable
    class(file_variable_t), intent(in) :: this

    units = this%units_

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the MUSICA units for the variable
  !!
  !! Unspecified file variable units are assumed to be in MUSICA units
  !!
  subroutine set_musica_units( this, units )

    use musica_assert,                 only : assert_msg

    !> File variable
    class(file_variable_t), intent(inout) :: this
    !> MUSICA units
    character(len=*), intent(in) :: units

    type(string_t) :: l_units

    if( this%musica_units_ .eq. "" ) this%musica_units_ = units
    l_units = units
    call assert_msg( 208760654, this%musica_units_%to_lower( ) .eq.           &
                                l_units%to_lower( ),                          &
                     "Attempting to change units for variable '"//            &
                     this%name_%to_char( )//"' from '"//                      &
                     this%musica_units_%to_char( )//"' to '"//                &
                     l_units%to_char( )//"'." )
    if( this%units_ .eq. "" ) this%units_ = this%musica_units_
    this%converter_ = convert_t( this%musica_units_, this%units_ )

  end subroutine set_musica_units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale, offset, convert units, and shift results for a set of data
  subroutine convert_to_musica_values( this, values )

    !> File variable
    class(file_variable_t), intent(in) :: this
    !> Data to process
    real(kind=musica_dk), intent(inout) :: values(:)

    integer(kind=musica_ik) :: i_val

    do i_val = 1, size( values )
      values( i_val ) = values( i_val ) * this%scale_factor_ + this%offset_
      values( i_val ) = this%converter_%to_standard( values( i_val ) )
      values( i_val ) = values( i_val ) + this%shift_
    end do

  end subroutine convert_to_musica_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale, offset, convert units, and shift a set of data to file values
  !! (Should only be called by extending types)
  subroutine convert_to_file_values( this, values )

    !> File variable
    class(file_variable_t), intent(in) :: this
    !> Data to process
    real(kind=musica_dk), intent(inout) :: values(:)

    integer(kind=musica_ik) :: i_val

    do i_val = 1, size( values )
      values( i_val ) = values( i_val ) - this%shift_
      values( i_val ) = this%converter_%to_non_standard( values( i_val ) )
      values( i_val ) = ( values( i_val ) - this%offset_ ) / this%scale_factor_
    end do

  end subroutine convert_to_file_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the properties of the variable
  subroutine do_print( this )

    !> File variable
    class(file_variable_t), intent(in) :: this

    write(*,*) "*** Variable: "//this%name_%to_char( )//" ***"
    write(*,*) "MUSICA name: "//this%musica_name_%to_char( )
    write(*,*) "units: "//this%units_%to_char( )
    write(*,*) "scale factor:", this%scale_factor_
    write(*,*) "offset:", this%offset_
    write(*,*) "shift:", this%shift_

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Equality comparison
  logical elemental function equals( a, b )

    !> File variable a
    class(file_variable_t), intent(in) :: a
    !> File variable b
    class(file_variable_t), intent(in) :: b

    equals = a%is_same_as( b )

  end function equals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Inequality comparison
  logical elemental function not_equals( a, b )

    !> File variable a
    class(file_variable_t), intent(in) :: a
    !> File variable b
    class(file_variable_t), intent(in) :: b

    not_equals = .not. a%is_same_as( b )

  end function not_equals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Loads matching and conversion options specified in the configuration data
  subroutine load_configured_options( this, file, config )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_file,                   only : file_t

    !> File variable
    class(file_variable_t), intent(inout) :: this
    !> Input/Output File
    class(file_t), intent(inout) :: file
    !> Configuration describing how to match to MUSICA variables
    !!
    !! If omitted, standard matching is applied
    type(config_t), intent(inout), optional :: config

    character(len=*), parameter :: my_name = "File variable matching"
    type(config_t) :: vars, var_data, shift_data
    type(string_t) :: default_variability, variability
    logical :: found, general_match

    ! default to File variable name
    this%musica_name_ = this%name_

    ! get the default variability (input variables only)
    if( file%is_input( ) ) then
      call config%get( "default variability", default_variability, my_name )
    end if

    ! get specific property matching if present
    call config%get( "properties", vars, my_name, found = found )

    ! look for specific and then general variable information
    general_match = .false.
    if( found ) then
      call vars%get( this%name_%to_char( ), var_data, my_name, found = found )
      if( .not. found ) then
        call vars%get( "*", var_data, my_name, found = found )
        general_match = found
      end if
    end if

    ! update matching criteria as specified in configuration
    if( found ) then
      call var_data%get( "MusicBox name", this%musica_name_, my_name,         &
                         default = this%musica_name_%to_char( ) )
      call var_data%get( "units", this%units_, my_name,                       &
                         default = this%units_%to_char( ) )
      if( file%is_input( ) ) then
        call var_data%get( "variability", variability, my_name,               &
                           default = default_variability%to_char( ) )
        if( variability .eq. "prognosed" ) then
          this%variability_ = VARIABILITY_PROGNOSED
        else if( variability .eq. "tethered" ) then
          this%variability_ = VARIABILITY_TETHERED
        else if( variability .eq. "fixed" ) then
          this%variability_ = VARIABILITY_FIXED
        else
          call die_msg( 180506671, "Invalid variability specified for "//     &
                        "file variable '"//this%name_%to_char( )//"': '"//    &
                        variability%to_char( )//"'" )
        end if
      end if
      call var_data%get( "shift first entry to", shift_data, my_name,         &
                         found = found )
      if( found ) then
        call this%set_shift( file, shift_data )
      end if
      if( general_match ) then
        this%musica_name_ = this%musica_name_%replace( "*",                   &
                                                       this%name_%to_char( ) )
      end if
    end if

    ! do standard name conversions
    call this%do_standard_name_conversions( )

  end subroutine load_configured_options

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does standard name conversions between file and MUSICA variables
  subroutine do_standard_name_conversions( this )

    !> File variable
    class(file_variable_t), intent(inout) :: this

    type(string_t), allocatable :: strs(:)

    associate( str => this%musica_name_ )
      str = str%replace( "CONC.", "chemical_species%" )
      str = str%replace( "ENV.",  "" )
      str = str%replace( "EMIS.", "emission_rates%" )
      str = str%replace( "LOSS.", "loss_rate_constants%" )
      str = str%replace( "PHOT.", "photolysis_rate_constants%" )
    end associate

    ! extract units if included in name
    strs = this%musica_name_%split( "." )
    if( size( strs ) .eq. 2 ) then
      this%musica_name_ = strs(1)
      this%units_       = strs(2)
    end if

  end subroutine do_standard_name_conversions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets a shift in input/output data
  !!
  !! (Currently only time shifts are supported.)
  subroutine set_shift( this, file, config )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_datetime,               only : datetime_t
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t

    !> File variable
    class(file_variable_t), intent(inout) :: this
    !> Input/Output File
    class(file_t), intent(inout) :: file
    !> Configuration data
    type(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name = "File variable shift"
    real(kind=musica_dk) :: values(1)
    type(datetime_t) :: shift
    type(file_dimension_range_t), allocatable :: var_dims(:)

    ! for now, just handle time shifts
    call this%set_musica_units( "s" )
    var_dims = this%get_dimensions( )
    call assert_msg( 802180458, size( var_dims ) .eq. 1,                      &
                     "Cannot shift multi-dimensional variable '"//            &
                     this%name_%to_char( )//"'" )
    call var_dims(1)%set( 1, 1 )
    call this%get_data( file, var_dims(1:1), values )
    shift = datetime_t( config )
    this%shift_ = shift%in_seconds( ) - values(1)

  end subroutine set_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Private constructor for common data elements
  !! (Should only be called by constructors of extending types)
  subroutine private_constructor( this, config, file, name, dimension_names,  &
      units )

    use musica_config,                 only : config_t
    use musica_file,                   only : file_t

    !> File variable
    class(file_variable_t), intent(inout) :: this
    !> Variable configuration
    type(config_t), intent(inout) :: config
    !> Input/Output File
    class(file_t), intent(inout) :: file
    !> Name used in the file for the variable
    type(string_t), intent(in) :: name
    !> Dimension names
    type(string_t), intent(in) :: dimension_names(:)
    !> Units used in the file for the variable
    !! (If not included, standard MUSICA units will be assumed)
    type(string_t), intent(in), optional :: units

    type(string_t) :: std_units

    this%name_            = name
    this%dimension_names_ = dimension_names
    if( present( units ) ) then
      this%units_ = units
    else
      this%units_ = ""
    end if
    this%musica_units_ = ""
    call this%load_configured_options( file, config )

  end subroutine private_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finds a File variable by name in a set of variables
  !!
  !! Variable matching is case-insensitive
  !!
  function find_variable_by_name( set, name, found ) result( var_id )

    use musica_assert,                 only : die

    !> Index of variable in set
    integer(musica_ik) :: var_id
    !> Set of File variables
    class(file_variable_ptr), intent(in) :: set(:)
    !> Variable name to locate
    type(string_t), intent(in) :: name
    !> Flag indicating whether variable was found
    logical, intent(out), optional :: found

    type(string_t) :: l_name, var_name

    l_name = name%to_lower( )
    do var_id = 1, size( set )
      if( set( var_id )%val_%name_%to_lower( ) .eq. l_name ) then
        if( present( found ) ) found = .true.
        return
      end if
    end do
    if( present( found ) ) then
      found = .false.
    else
      call die( 249940482 )
    end if

  end function find_variable_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finds a File variable id by its expected MUSICA name
  !!
  !! Variable matching is case-insensitive
  !!
  function find_variable_by_musica_name( set, name, found ) result( var_id )

    use musica_assert,                 only : die

    !> Index of variable in set
    integer(musica_ik) :: var_id
    !> Set of File variables
    class(file_variable_ptr), intent(in) :: set(:)
    !> Domain variable name to locate
    type(string_t), intent(in) :: name
    !> Flag indicating whether variable was found
    logical, intent(out), optional :: found

    type(string_t) :: l_name, musica_name

    l_name = name%to_lower( )
    do var_id = 1, size( set )
      if( set( var_id )%val_%musica_name_%to_lower( ) .eq. l_name ) then
        if( present( found ) ) found = .true.
        return
      end if
    end do
    if( present( found ) ) then
      found = .false.
    else
      call die( 464048558 )
    end if

  end function find_variable_by_musica_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a unique file variable pointer
  elemental subroutine file_variable_ptr_finalize( this )

    !> File variable pointer
    type(file_variable_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine file_variable_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_variable
