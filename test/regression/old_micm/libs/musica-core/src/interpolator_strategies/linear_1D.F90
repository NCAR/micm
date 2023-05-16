! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_interpolator_linear module

!> Linear strategy for 1D interpolators
module musica_interpolator_linear_1D

  implicit none
  private

  public :: strategy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Perform a linear 1D interpolation
  !!
  !! Note: not optimized. Use during initialization or consider optimizing.
  !!
  !! An example of how to use the interpolator object, which uses the linear
  !! 1D interpolation strategy follows:
  !!
  !! \snippet test/interpolator_strategies/linear_1D.F90 Linear interpolator example
  !!
  !! Output:
  !! \code{bash}
  !!    5.0000000000000000        5.0000000000000000        10.000000000000000        10.000000000000000        15.000000000000000        15.000000000000000
  !! \endcode
  function strategy( from_grid, to_grid ) result( map )

    use musica_constants,              only : musica_dk
    use musica_grid,                   only : grid_t, grid_iterator_t,        &
                                              grid_section_t
    use musica_iterator,               only : iterator_t
    use musica_interpolator,           only : interpolator_element_t

    class(interpolator_element_t), allocatable :: map(:)
    class(grid_t), intent(in) :: from_grid
    class(grid_t), intent(in) :: to_grid

    class(grid_iterator_t), pointer :: iter_from, iter_to
    class(grid_section_t), allocatable :: from_section, to_section
    type(interpolator_element_t), allocatable :: local_map(:)
    type(interpolator_element_t) :: element
    real(kind=musica_dk) :: overlap

    allocate( local_map(0) )
    iter_from => from_grid%iterator( )
    iter_to   => to_grid%iterator( )
    do while( iter_to%next( ) )
      to_section   = to_grid%section( iter_to )
      do while( iter_from%next( ) )
        from_section = from_grid%section( iter_from )
        overlap = from_section%overlap( to_section )
        if( overlap .gt. 0.0 ) then
          element = interpolator_element_t( iter_from, iter_to,               &
                                            overlap / from_section%range( ) )
          local_map = [ local_map, element ]
        end if
      end do
      call iter_from%reset( )
    end do
    allocate( map, source = local_map )
    deallocate( iter_from )
    deallocate( iter_to   )

  end function strategy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_interpolator_linear_1D
