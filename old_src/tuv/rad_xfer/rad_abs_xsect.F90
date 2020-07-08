module rad_abs_xsect

  use phot_kind_mod, only: rk => kind_phot

  implicit none

  public

  real(rk), protected, allocatable :: o2_xs(:)
  real(rk), protected, allocatable :: so2_xs(:)
  real(rk), protected, allocatable :: wl(:)
  real(rk), protected, allocatable :: wc(:)

  integer, protected :: nwave

contains

  subroutine rad_abs_xsect_init( filepath, errmsg, errflg )

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
       errmsg = 'rad_abs_xsect_init: failed to open '//trim(filepath)
       return
    end if

    ! get dimensions
    ret = nf90_inq_dimid( ncid, 'nwave', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get nwave id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=nwave )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get nwave'
       return
    end if

    ! allocate memory
    allocate( wl(nwave+1), wc(nwave), o2_xs(nwave), so2_xs(nwave), stat=astat )
    if( astat /= 0 ) then
       errflg = astat
       errmsg = 'rad_abs_xsect_init: failed to allocate memory'
       return
    end if

    ! read xsect data
    ret = nf90_inq_varid( ncid, 'o2_xs', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get o2_xs variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, o2_xs )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to read o2_xs variable'
       return
    end if
    ret = nf90_inq_varid( ncid, 'so2_xs', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get so2_xs variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, so2_xs )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to read so2_xs variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'wl', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get wl variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, wl )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to read wl variable'
       return
    end if
    ret = nf90_inq_varid( ncid, 'wc', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get wc variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, wc )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to read wc variable'
       return
    end if

    ! close the file
    ret = nf90_close( ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to close '//trim(filepath)
       return
    end if

  end subroutine rad_abs_xsect_init

end module rad_abs_xsect
