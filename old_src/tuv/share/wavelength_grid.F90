module wavelength_grid

  use phot_kind_mod, only: rk => kind_phot

  implicit none

  public

  real(rk), protected, allocatable :: wl(:)
  real(rk), protected, allocatable :: wc(:)
  integer, protected :: nwave
  
contains

  subroutine wavelength_grid_init( filepath, errmsg, errflg )

    use netcdf

    character(len=*), intent(in ) :: filepath
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: ncid, dimid, varid
    integer :: astat, ret

    errmsg = ' '
    errflg = 0

    ! open file
    ret = nf90_open( trim(filepath), nf90_noclobber, ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'wavelength_grid_init: failed to open '//trim(filepath)
       return
    end if

    ! get dimensions
    ret = nf90_inq_dimid( ncid, 'nwave', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'wavelength_grid_init: failed to get nwave id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=nwave )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'wavelength_grid_init: failed to get nwave'
       return
    end if

    ! allocate memory
    allocate( wl(nwave+1), wc(nwave), stat=astat )
    if( astat /= 0 ) then
       errflg = astat
       errmsg = 'wavelength_grid_init: failed to allocate memory'
       return
    end if

    ret = nf90_inq_varid( ncid, 'wl', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'wavelength_grid_init: failed to get wl variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, wl )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'wavelength_grid_init: failed to read wl variable'
       return
    end if
    ret = nf90_inq_varid( ncid, 'wc', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'wavelength_grid_init: failed to get wc variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, wc )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'wavelength_grid_init: failed to read wc variable'
       return
    end if

    ! close the file
    ret = nf90_close( ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'wavelength_grid_init: failed to close '//trim(filepath)
       return
    end if

  end subroutine wavelength_grid_init

end module wavelength_grid
