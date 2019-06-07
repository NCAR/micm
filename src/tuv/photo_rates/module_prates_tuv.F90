module  module_prates_tuv
  use phot_kind_mod, only: rk => kind_phot
  use module_rxn, only : xsqy_table => xsqy_tab, the_subs, npht_tab, rxn_init
  use module_rxn, only : get_initialization, set_initialization

  implicit none

  private
  public :: calc_tuv_init
  public :: calc_tuv_prates
  public :: rxn_ndx
  public :: get_xsqy_tab
  public :: nwave

  integer, protected :: nwave
  integer :: nconc, ntemp

  real(rk), allocatable :: z_temp_data(:), temp_data(:)
  real(rk), allocatable :: xsqy_tab(:,:,:,:)
  real(rk), allocatable :: z_o3_data(:), o3_data(:)
  real(rk), allocatable :: z_air_dens_data(:), air_dens_data(:)
  real(rk), allocatable :: wl(:)
  real(rk), allocatable :: wc(:)
  real(rk), allocatable :: dw(:)
  real(rk), allocatable :: w_fac(:)
  real(rk), allocatable :: etfl(:)
  real(rk), allocatable :: temp_tab(:)
  real(rk), allocatable :: conc_tab(:)
  real(rk), allocatable :: del_temp_tab(:)
  real(rk), allocatable :: del_conc_tab(:)
  real(rk), allocatable :: o3_xs_tab(:,:)
  real(rk), allocatable :: no2_xs_tab(:,:)
  character(len=32), allocatable :: tuv_jname(:)
  integer :: n_temp_data, n_o3_data, n_air_dens_data

  integer :: nj 
  integer :: j_o2_ndx = 1
  logical :: is_full_tuv = .true.

  logical, allocatable :: xsqy_is_zdep(:)
  integer, protected, allocatable :: rxn_ndx(:)

  integer, parameter :: nlambda_start = 1
  real(rk)    :: esfact = 1.0_rk
contains

  !------------------------------------------------------------------------------
  ! initialize
  !------------------------------------------------------------------------------
  subroutine calc_tuv_init( full_tuv, jnames, errmsg, errflg )

    logical, intent(in) :: full_tuv
    character(len=*), intent(in) :: jnames(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    
    integer :: i,j,n

    is_full_tuv = full_tuv
    nj = size(jnames)
    allocate( tuv_jname(nj) )
    tuv_jname(:) = jnames(:)
    
    if( .not. is_full_tuv ) then
       allocate( xsqy_is_zdep(nj) )
       xsqy_is_zdep(:) = .false.
       if( j_o2_ndx > 0 ) then
          xsqy_is_zdep(j_o2_ndx) = .true.
       endif
       do n = 1,nj
          if( n /= j_o2_ndx ) then
             t_loop :      do j = 1,nconc
                do i = 1,ntemp
                   if( any( xsqy_tab(:,i,j,n) /= xsqy_tab(:,1,1,n) ) ) then
                      xsqy_is_zdep(n) = .true.
                      exit t_loop
                   endif
                end do
             end do t_loop
          endif
       end do
    else
       call rxn_init( nwave+1,wl, errmsg, errflg )
       allocate( rxn_ndx(nj) )
       rxn_ndx(1:nj) = -1
       do j = 1,nj
          if( j /= j_o2_ndx ) then
             do n = 2,npht_tab
                if( trim(xsqy_table(n)%wrf_label) == trim(tuv_jname(j)) ) then
                   rxn_ndx(j) = n
                   exit
                endif
             enddo
          endif
       enddo
    endif       

  end subroutine calc_tuv_init

  !------------------------------------------------------------------------------
  ! compute Js for a column given photon fluxes in each grid box of the column (nlevs)
  !------------------------------------------------------------------------------
  subroutine calc_tuv_prates(kts,kte,nlevs, tlev, dens_air, rad_fld, srb_o2_xs, tuv_prate, errmsg, errflg)

    USE,INTRINSIC :: IEEE_ARITHMETIC

    ! Args
    integer, intent(in) :: kts,kte,nlevs
    real(rk), intent(in)    :: tlev(kts:kte)     ! K
    real(rk), intent(in)    :: dens_air(kts:kte) ! molecules / cm3
    real(rk), intent(in)    :: rad_fld(nwave,kts:kte)
    real(rk), intent(in)    :: srb_o2_xs(nwave,kts:kte)
    real(rk), intent(out)   :: tuv_prate(nlevs, nj) ! /sec

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Locals
    real(rk) :: xsect(nwave)
    real(rk) :: xsqy(nwave,kts:kte)
    real(rk) :: rad_fld_tpose(kts:kte,nwave)

    integer :: n,jndx
    integer :: k
    real(rk) :: dummy(nlevs)
    logical :: rxn_initialized
    real(rk) :: xnan

    xnan =  IEEE_VALUE(xnan,IEEE_QUIET_NAN)
    tuv_prate = xnan
    rad_fld_tpose = xnan

    if( .not. is_full_tuv ) then
       if( any( .not. xsqy_is_zdep(:) ) ) then
          rad_fld_tpose = transpose( rad_fld )
       endif
    elseif( any( xsqy_table(1:nj)%tpflag == 0 ) ) then
       rad_fld_tpose = transpose( rad_fld )
    endif

    if( is_full_tuv ) then
       rxn_initialized = .not. get_initialization()
       if( .not. rxn_initialized ) then
          do n = 1,nj
             jndx = rxn_ndx(n)
             if( jndx /= -1 ) then
                call the_subs(jndx)%xsqy_sub(nwave+1,wl,wc,nlevs,dummy,dummy,jndx, errmsg, errflg )
             endif
          enddo
          call set_initialization( status=.false. )
       endif
    endif

    rate_loop: do n = 1,nj
    xsect = xnan
    xsqy = xnan
       !---------------------------------------------------------------------
       ! set cross-section x quantum yields
       !---------------------------------------------------------------------
       if( n /= j_o2_ndx ) then
          if( .not. is_full_tuv ) then
             if( xsqy_is_zdep(n) ) then
                call xsqy_int( n, xsqy, tlev(kts:kte), dens_air(kts:kte) )
             endif
          else
             jndx = rxn_ndx(n)
             if( jndx /= -1 ) then
                if( xsqy_table(jndx)%tpflag /= 0 ) then
                   call the_subs(jndx)%xsqy_sub(nwave+1,wl,wc,nlevs,tlev,dens_air,jndx, errmsg, errflg)
                endif
             endif
          endif
       elseif( .not. is_full_tuv ) then
          call sjo2( kte, nwave, srb_o2_xs, xsqy )
       endif
       !---------------------------------------------------------------------
       ! compute tuv photorates
       !---------------------------------------------------------------------
       if( .not. is_full_tuv ) then
          if( xsqy_is_zdep(n) ) then
             do k = kts,kte
                xsect(nlambda_start:nwave) = xsqy(nlambda_start:nwave,k)*w_fac(nlambda_start:nwave)*esfact
                tuv_prate(k,n) = dot_product( rad_fld(nlambda_start:nwave,k),xsect(nlambda_start:nwave) )
             end do
          else
             xsect(nlambda_start:nwave) = xsqy_tab(nlambda_start:nwave,1,1,n)*w_fac(nlambda_start:nwave)*esfact
             tuv_prate(:,n) = matmul( rad_fld_tpose(:,nlambda_start:nwave),xsect(nlambda_start:nwave) )
          endif
       else
          if( n /= j_o2_ndx ) then
             if( xsqy_table(jndx)%tpflag > 0 ) then
                do k = kts,kte
                   xsect(nlambda_start:nwave) = xsqy_table(jndx)%sq(k,nlambda_start:nwave)*w_fac(nlambda_start:nwave)*esfact
                   tuv_prate(k,n) = dot_product( rad_fld(nlambda_start:nwave,k),xsect(nlambda_start:nwave) )
                end do
             else
                xsect(nlambda_start:nwave) = xsqy_table(jndx)%sq(nlambda_start:nwave,1)*w_fac(nlambda_start:nwave)*esfact
                tuv_prate(:,n) = matmul( rad_fld_tpose(:,nlambda_start:nwave),xsect(nlambda_start:nwave) )
             endif
          else
             do k = kts,kte
                xsect(nlambda_start:nwave) = srb_o2_xs(nlambda_start:nwave,k)*w_fac(nlambda_start:nwave)*esfact
                tuv_prate(k,n) = dot_product( rad_fld(nlambda_start:nwave,k),xsect(nlambda_start:nwave) )
             end do
          endif
       endif
    end do rate_loop

  end subroutine calc_tuv_prates

  subroutine get_xsqy_tab(xsqy_filepath, errmsg, errflg )
    !---------------------------------------------------------------------
    !	... read in the cross section,quantum yield tables
    !---------------------------------------------------------------------

    use params_mod, only : hc
    use netcdf

    character(len=*), intent(in)  :: xsqy_filepath
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !---------------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------------
    integer :: astat, ierr, ret
    integer :: m
    integer :: ncid, dimid, varid
    character(len=64) :: varname

    ret = nf90_open( trim(xsqy_filepath), nf90_noclobber, ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to open file ' // trim(xsqy_filepath)
       return
    end if

    !---------------------------------------------------------------------
    !	... get dimensions
    !---------------------------------------------------------------------
    ret = nf90_inq_dimid( ncid, 'nwave', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get nwave id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=nwave )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get nwave'
       return
    end if
    ret = nf90_inq_dimid( ncid, 'ntemp', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get ntemp id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=ntemp )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get ntemp'
       return
    end if
    ret = nf90_inq_dimid( ncid, 'nconc', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get nconc id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=nconc )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get nconc'
       return
    end if
    ret = nf90_inq_dimid( ncid, 'n_temp_data', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get n_temp_data id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=n_temp_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get n_temp_data'
       return
    end if
    ret = nf90_inq_dimid( ncid, 'n_o3_data', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get n_o3_data id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=n_o3_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get n_o3_data'
       return
    end if
    ret = nf90_inq_dimid( ncid, 'n_air_dens_data', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get n_air_dens_data id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=n_air_dens_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get n_air_dens_data'
       return
    end if

    !---------------------------------------------------------------------
    !	... allocate module arrays
    !---------------------------------------------------------------------
    ierr = 0
    allocate( z_temp_data(n_temp_data), z_o3_data(n_o3_data), &
         z_air_dens_data(n_air_dens_data),stat=astat )
    ierr = astat + ierr
    allocate( temp_data(n_temp_data), o3_data(n_o3_data), &
         air_dens_data(n_air_dens_data),stat=astat )
    ierr = astat + ierr
    allocate( wl(nwave+1), wc(nwave), dw(nwave), w_fac(nwave), &
         etfl(nwave), stat=astat )
    ierr = astat + ierr
    if( .not. is_full_tuv ) then
       allocate( temp_tab(ntemp), conc_tab(nconc), stat=astat )
       ierr = astat + ierr
       allocate( del_temp_tab(ntemp-1), del_conc_tab(nconc-1), stat=astat )
       ierr = astat + ierr
    endif
    allocate( o3_xs_tab(nwave,ntemp), no2_xs_tab(nwave,ntemp), stat=astat )
    ierr = astat + ierr
    if( .not. is_full_tuv ) then
       allocate( xsqy_tab(nwave,ntemp,nconc,nj), stat=astat )
       ierr = astat + ierr
    endif
    if( ierr /= 0 ) then
       errmsg = 'get_xsqy_tab: failed to allocate'
       errflg = ierr
       return
    endif
    
    !---------------------------------------------------------------------
    !	... read arrays
    !---------------------------------------------------------------------
    ret = nf90_inq_varid( ncid, 'z_temp_data', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get z_temp_data variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, z_temp_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read z_temp_data variable'
       return
    end if


    ret = nf90_inq_varid( ncid, 'z_o3_data', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get z_o3_data variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, z_o3_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read z_o3_data variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'z_air_dens_data', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get z_air_dens_data variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, z_air_dens_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read z_air_dens_data variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'temp_data', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get temp_data variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, temp_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read temp_data variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'o3_data', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get o3_data variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, o3_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read o3_data variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'air_dens_data', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get air_dens_data variable id'
       return
    end if


    ret = nf90_get_var( ncid, varid, air_dens_data )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read air_dens_data variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'wl', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get wl variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, wl )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read wl variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'wc', varid )
     if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get wc variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, wc )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read wc variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'etf', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get etfl variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, etfl )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read etfl variable'
       return
    end if

    if( .not. is_full_tuv ) then
       ret = nf90_inq_varid( ncid, 'temps', varid )
       if( ret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'get_xsqy_tab: failed to temp_tab variable id'
          return
       end if

       ret = nf90_get_var( ncid, varid, temp_tab )
       if( ret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'get_xsqy_tab: failed to read temp_tab variable'
          return
       end if

       ret = nf90_inq_varid( ncid, 'concs', varid )
       if( ret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'get_xsqy_tab: failed to conc_tab variable id'
          return
       end if

       ret = nf90_get_var( ncid, varid, conc_tab )
       if( ret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'get_xsqy_tab: failed to read conc_tab variable'
          return
       end if

    endif
    
    ret = nf90_inq_varid( ncid, 'o3_xs', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to o3_xs_tab variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, o3_xs_tab )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read o3_xs_tab variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'no2_xs', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to no2_xs_tab variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, no2_xs_tab )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read no2_xs_tab variable'
       return
    end if
    if( .not. is_full_tuv ) then
       do m = 1,nj
          varname = trim(tuv_jname(m)) // '_xsqy'
          ret = nf90_inq_varid( ncid, trim(varname), varid )
          if( ret /= nf90_noerr ) then
             errflg = 1
             errmsg = 'get_xsqy_tab: failed to ' // trim(varname) //' variable id'
             return
          end if
          ret = nf90_get_var( ncid, varid, xsqy_tab(:,:,:,m) )
          if( ret /= nf90_noerr ) then
             errflg = 1
             errmsg = 'get_xsqy_tab: failed to read ' // trim(varname) // ' variable'
             return
          end if
       end do
    endif

    if( .not. is_full_tuv ) then
       del_temp_tab(:ntemp-1) = 1._rk/(temp_tab(2:ntemp) - temp_tab(1:ntemp-1))
       del_conc_tab(:nconc-1) = 1._rk/(conc_tab(2:nconc) - conc_tab(1:nconc-1))
    endif
    dw(:nwave)    = wl(2:nwave+1) - wl(1:nwave)
    w_fac(:nwave) = dw(:nwave)*etfl(:nwave)*1.e-13_rk*wc(:nwave)/hc

    ret = nf90_close( ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to close file ' // trim(xsqy_filepath)
       return
    end if

  end subroutine get_xsqy_tab

      subroutine xsqy_int( n, xsqy, tlev, dens_air )
!---------------------------------------------------------------------
!	... interpolate m,t tables for xs * qy
!---------------------------------------------------------------------

      integer, intent(in)  :: n 
      real(rk),    intent(in)  :: tlev(:)
      real(rk),    intent(in)  :: dens_air(:)
      real(rk),    intent(out) :: xsqy(:,:)

      real(rk), parameter :: m0 = 2.45e19_rk
      integer :: tndx, mndx, tndxp1, mndxp1
      integer :: k, ku
      real(rk)    :: temp, dens
      real(rk)    :: w(4)
      real(rk)    :: del_t, del_d

      ku = size( tlev )
      do k = 1,ku
        temp = tlev(k)
        do tndx = 1,ntemp
          if( temp_tab(tndx) > temp ) then
            exit
          endif
        end do
        tndx = max( min( tndx,ntemp ) - 1,1 )
        tndxp1 = tndx + 1
        del_t = max( 0._rk,min( 1._rk,(temp - temp_tab(tndx))*del_temp_tab(tndx) ) )

!       dens = dens_air(k)
        dens = dens_air(k)/m0
        do mndx = 1,nconc
          if( conc_tab(mndx) > dens ) then
            exit
          endif
        end do
        mndx = max( min( mndx,nconc ) - 1,1 )
        mndxp1 = mndx + 1
        del_d = max( 0._rk,min( 1._rk,(dens - conc_tab(mndx))*del_conc_tab(mndx) ) )

        w(1) = (1._rk - del_t)*(1._rk - del_d)
        w(2) = del_t*(1._rk - del_d)
        w(3) = (1._rk - del_t)*del_d
        w(4) = del_t*del_d

        xsqy(1:nwave,k) = w(1)*xsqy_tab(1:nwave,tndx,mndx,n) &
                        + w(2)*xsqy_tab(1:nwave,tndxp1,mndx,n) &
                        + w(3)*xsqy_tab(1:nwave,tndx,mndxp1,n) &
                        + w(4)*xsqy_tab(1:nwave,tndxp1,mndxp1,n)
      end do

      end subroutine xsqy_int

      subroutine xs_int( xs, tlev, xs_tab )
!---------------------------------------------------------------------
!	... interpolate tables for xs
!---------------------------------------------------------------------

      real(rk),    intent(in)  :: tlev(:)
      real(rk),    intent(in)  :: xs_tab(:,:)
      real(rk),    intent(out) :: xs(:,:)

      integer :: tndx, tndxp1
      integer :: k, ku
      real(rk)    :: temp
      real(rk)    :: w(2)
      real(rk)    :: del_t

      ku = size( tlev )
      do k = 1,ku
        temp = tlev(k)
        do tndx = 1,ntemp
          if( temp_tab(tndx) > temp ) then
            exit
          endif
        end do
        tndx = max( min( tndx,ntemp ) - 1,1 )
        tndxp1 = tndx + 1
        del_t = max( 0._rk,min( 1._rk,(temp - temp_tab(tndx))*del_temp_tab(tndx) ) )

        w(1) = (1._rk - del_t)
        w(2) = del_t

        xs(1:nwave,k) = w(1)*xs_tab(1:nwave,tndx) &
                      + w(2)*xs_tab(1:nwave,tndxp1)
      end do

      end subroutine xs_int

      SUBROUTINE sjo2( nlyr, nwave, xso2, xsqy )
!-----------------------------------------------------------------------------
!=  PURPOSE:
!=  Update the weighting function (cross section x quantum yield) for O2
!=  photolysis.  The strong spectral variations in the O2 cross sections are
!=  parameterized into a few bands for Lyman-alpha (121.4-121.9 nm, one band)
!=  and Schumann-Runge (174.4-205.8, nsrb bands) regions. The parameterizations
!=  depend on the overhead O2 column, and therefore on altitude and solar
!=  zenith angle, so they need to be updated at each time/zenith step.
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)
!=           wavelength grid
!=  XSO2   - REAL, molecular absorption cross section in SR bands at      (I)
!=           each specified altitude and wavelength.  Includes Herzberg
!=            continuum.
!=  NJ     - INTEGER, index of O2 photolysis in array SQ                  (I)
!=  xsqy   - REAL, cross section x quantum yield (cm^2) for each          (O)
!=           photolysis reaction, at each wavelength and each altitude level
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      INTEGER, intent(in)    :: nlyr, nwave
      REAL(rk),    intent(in)    :: xso2(:,:)
      REAL(rk),    intent(inout) :: xsqy(:,:)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      INTEGER :: k

!-----------------------------------------------------------------------------
! O2 + hv -> O + O
! quantum yield assumed to be unity
! assign cross section values at all wavelengths and at all altitudes
!      qy = 1.
!-----------------------------------------------------------------------------
      DO k = 1, nlyr
        xsqy(:nwave,k) = xso2(:nwave,k)
      END DO

      END SUBROUTINE sjo2


end module module_prates_tuv
