module  module_prates_tuv
  use module_rxn, only : xsqy_table => xsqy_tab, the_subs, npht_tab, rxn_init
  use module_rxn, only : get_initialization, set_initialization
  use params_mod, only : m2s
  use la_srb_mod, only : sjo2, init_srb

  implicit none

  integer :: nconc, ntemp, nwave

  real, allocatable :: z_temp_data(:), temp_data(:)
  real, allocatable :: xsqy_tab(:,:,:,:)
  real, allocatable :: z_o3_data(:), o3_data(:)
  real, allocatable :: z_air_dens_data(:), air_dens_data(:)
  real, allocatable :: wl(:)
  real, allocatable :: wc(:)
  real, allocatable :: dw(:)
  real, allocatable :: w_fac(:)
  real, allocatable :: etfl(:)
  real, allocatable :: temp_tab(:)
  real, allocatable :: conc_tab(:)
  real, allocatable :: del_temp_tab(:)
  real, allocatable :: del_conc_tab(:)
  real, allocatable :: o2_xs(:)
  real, allocatable :: so2_xs(:)
  real, allocatable :: o3_xs_tab(:,:)
  real, allocatable :: no2_xs_tab(:,:)
  character(len=32), allocatable :: tuv_jname(:)
  integer :: n_temp_data, n_o3_data, n_air_dens_data

  integer :: nj 
  integer :: j_o2_ndx = 1
  logical :: is_full_tuv = .true.

  logical, allocatable :: xsqy_is_zdep(:)
  integer, allocatable :: rxn_ndx(:)

  integer :: nlambda_start = 1
  real    :: esfact = 1.0
contains

  !------------------------------------------------------------------------------
  ! initialize
  !------------------------------------------------------------------------------
  subroutine calc_tuv_init( xsqy_filepath, full_tuv, jnames )

    character(len=*), intent(in) :: xsqy_filepath
    logical, intent(in) :: full_tuv
    character(len=*), intent(in) :: jnames(:)
    
    integer :: i,j,n

    is_full_tuv = full_tuv
    nj = size(jnames)
    allocate( tuv_jname(nj) )
    tuv_jname(:) = jnames(:)
    
    call get_xsqy_tab(xsqy_filepath)
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
       call rxn_init( nwave+1,wl )
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

    call init_srb()
  end subroutine calc_tuv_init

  !------------------------------------------------------------------------------
  ! compute Js for a column given photon fluxes in each grid box of the column (nlevs)
  !------------------------------------------------------------------------------
  subroutine calc_tuv_prates(kts,kte,nlevs, wl,wc, tlev,dens_air,rad_fld,srb_o2_xs, njs,tuv_prate)

    USE,INTRINSIC :: IEEE_ARITHMETIC

    ! Args
    integer, intent(in) :: kts,kte,nlevs
    real, intent(in)    :: wl(nwave+1), wc(nwave)
    real, intent(in)    :: tlev(kts:kte)     ! K
    real, intent(in)    :: dens_air(kts:kte) ! molecules / cm3
    real, intent(in)    :: rad_fld(nwave,kts:kte)
    real, intent(in)    :: srb_o2_xs(nwave,kts:kte)
    integer, intent(in) :: njs
    real, intent(out)   :: tuv_prate(nlevs, njs) ! /sec

    ! Locals
    real :: xsect(nwave)
    real :: xsqy(nwave,kts:kte)
    real :: rad_fld_tpose(kts:kte,nwave)

    integer :: n,jndx
    integer :: k
    real :: dummy(nlevs)
    logical :: rxn_initialized
    real :: xnan

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
                call the_subs(jndx)%xsqy_sub(nwave+1,wl,wc,nlevs,dummy,dummy,jndx)
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
                   call the_subs(jndx)%xsqy_sub(nwave+1,wl,wc,nlevs,tlev,dens_air,jndx)
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

  subroutine get_xsqy_tab(xsqy_filepath)
    !---------------------------------------------------------------------
    !	... read in the cross section,quantum yield tables
    !---------------------------------------------------------------------

    use params_mod, only : hc
    use la_srb_mod, only : ila, isrb
    use la_srb_mod, only : nchebev_term, nchebev_wave
    use la_srb_mod, only : chebev_ac, chebev_bc
    use netcdf
    use tuv_error_mod

    character(len=*), intent(in) :: xsqy_filepath

    !---------------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------------
    integer :: astat, ierr
    integer :: m
    integer :: ncid, dimid, varid
    character(len=132) :: err_msg
    character(len=64)  :: varname


    err_msg = 'get_xsqy_tab: failed to open file ' // trim(xsqy_filepath)
    call handle_ncerr( nf90_open( trim(xsqy_filepath), nf90_noclobber, ncid ), trim(err_msg) )
    !---------------------------------------------------------------------
    !	... get dimensions
    !---------------------------------------------------------------------
    err_msg = 'get_xsqy_tab: failed to get nwave id'
    call handle_ncerr( nf90_inq_dimid( ncid, 'nwave', dimid ), trim(err_msg) ) 
    err_msg = 'get_xsqy_tab: failed to get nwave'
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=nwave ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get ntemp id'
    call handle_ncerr( nf90_inq_dimid( ncid, 'ntemp', dimid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get ntemp'
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=ntemp ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get nconc id'
    call handle_ncerr( nf90_inq_dimid( ncid, 'nconc', dimid ), trim(err_msg) ) 
    err_msg = 'get_xsqy_tab: failed to get nconc'
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=nconc ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get nchebev_term id'
    call handle_ncerr( nf90_inq_dimid( ncid, 'nchebev_term', dimid ), trim(err_msg) ) 
    err_msg = 'get_xsqy_tab: failed to get nchebev'
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=nchebev_term ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get nchebev_wave id'
    call handle_ncerr( nf90_inq_dimid( ncid, 'nchebev_wave', dimid ), trim(err_msg) ) 
    err_msg = 'get_xsqy_tab: failed to get nchebev'
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=nchebev_wave ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get n_temp_data id'
    call handle_ncerr( nf90_inq_dimid( ncid, 'n_temp_data', dimid ), trim(err_msg) ) 
    err_msg = 'get_xsqy_tab: failed to get n_temp_data'
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=n_temp_data ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get n_o3_data id'
    call handle_ncerr( nf90_inq_dimid( ncid, 'n_o3_data', dimid ), trim(err_msg) ) 
    err_msg = 'get_xsqy_tab: failed to get n_o3_data'
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=n_o3_data ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get n_air_dens_data id'
    call handle_ncerr( nf90_inq_dimid( ncid, 'n_air_dens_data', dimid ), trim(err_msg) ) 
    err_msg = 'get_xsqy_tab: failed to get n_air_dens_data'
    call handle_ncerr( nf90_inquire_dimension( ncid, dimid, len=n_air_dens_data ), trim(err_msg) )

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
    allocate( chebev_ac(nchebev_term,nchebev_wave), stat=astat )
    ierr = astat + ierr
    allocate( chebev_bc(nchebev_term,nchebev_wave), stat=astat )
    ierr = astat + ierr
    allocate( o2_xs(nwave), so2_xs(nwave), stat=astat )
    ierr = astat + ierr
    allocate( o3_xs_tab(nwave,ntemp), no2_xs_tab(nwave,ntemp), stat=astat )
    ierr = astat + ierr
    if( .not. is_full_tuv ) then
       allocate( xsqy_tab(nwave,ntemp,nconc,nj), stat=astat )
       ierr = astat + ierr
    endif
    if( ierr /= 0 ) then
       call tuv_error_fatal( 'tuv_init: failed to allocate z_temp_data ...  xsqy_tab' )
    endif
    !---------------------------------------------------------------------
    !	... read arrays
    !---------------------------------------------------------------------
    err_msg = 'get_xsqy_tab: failed to get z_temp_data variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'z_temp_data', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read z_temp_data variable'
    call handle_ncerr( nf90_get_var( ncid, varid, z_temp_data ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get z_o3_data variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'z_o3_data', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read z_o3_data variable'
    call handle_ncerr( nf90_get_var( ncid, varid, z_o3_data ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get z_air_dens_data variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'z_air_dens_data', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read z_air_dens_data variable'
    call handle_ncerr( nf90_get_var( ncid, varid, z_air_dens_data ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get temp_data variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'temp_data', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read temp_data variable'
    call handle_ncerr( nf90_get_var( ncid, varid, temp_data ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get o3_data variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'o3_data', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read o3_data variable'
    call handle_ncerr( nf90_get_var( ncid, varid, o3_data ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get air_dens_data variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'air_dens_data', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read air_dens_data variable'
    call handle_ncerr( nf90_get_var( ncid, varid, air_dens_data ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get wl variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'wl', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read wl variable'
    call handle_ncerr( nf90_get_var( ncid, varid, wl ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get wc variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'wc', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read wc variable'
    call handle_ncerr( nf90_get_var( ncid, varid, wc ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get etfl variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'etf', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read etfl variable'
    call handle_ncerr( nf90_get_var( ncid, varid, etfl ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get chebev_ac variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'chebev_ac', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read chebev_ac variable'
    call handle_ncerr( nf90_get_var( ncid, varid, chebev_ac ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to get chebev_bc variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'chebev_bc', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read chebev_bc variable'
    call handle_ncerr( nf90_get_var( ncid, varid, chebev_bc ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get ila variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'ila', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read ila variable'
    call handle_ncerr( nf90_get_var( ncid, varid, ila ), trim(err_msg) )

    err_msg = 'get_xsqy_tab: failed to get isrb variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'isrb', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read isrb variable'
    call handle_ncerr( nf90_get_var( ncid, varid, isrb ), trim(err_msg) )

    if( .not. is_full_tuv ) then
       err_msg = 'get_xsqy_tab: failed to temp_tab variable id'
       call handle_ncerr( nf90_inq_varid( ncid, 'temps', varid ), trim(err_msg) )
       err_msg = 'get_xsqy_tab: failed to read temp_tab variable'
       call handle_ncerr( nf90_get_var( ncid, varid, temp_tab ), trim(err_msg) )
       err_msg = 'get_xsqy_tab: failed to conc_tab variable id'
       call handle_ncerr( nf90_inq_varid( ncid, 'concs', varid ), trim(err_msg) )
       err_msg = 'get_xsqy_tab: failed to read conc_tab variable'
       call handle_ncerr( nf90_get_var( ncid, varid, conc_tab ), trim(err_msg) )
    endif
    err_msg = 'get_xsqy_tab: failed to o2_xs variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'o2_xs', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read o2_xs variable'
    call handle_ncerr( nf90_get_var( ncid, varid, o2_xs ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to so2_xs variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'so2_xs', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read so2_xs variable'
    call handle_ncerr( nf90_get_var( ncid, varid, so2_xs ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to o3_xs_tab variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'o3_xs', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read o3_xs_tab variable'
    call handle_ncerr( nf90_get_var( ncid, varid, o3_xs_tab ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to no2_xs_tab variable id'
    call handle_ncerr( nf90_inq_varid( ncid, 'no2_xs', varid ), trim(err_msg) )
    err_msg = 'get_xsqy_tab: failed to read no2_xs_tab variable'
    call handle_ncerr( nf90_get_var( ncid, varid, no2_xs_tab ), trim(err_msg) )
    if( .not. is_full_tuv ) then
       do m = 1,nj
          varname = trim(tuv_jname(m)) // '_xsqy'
          err_msg = 'get_xsqy_tab: failed to ' // trim(varname) //' variable id'
          call handle_ncerr( nf90_inq_varid( ncid, trim(varname), varid ), trim(err_msg) )
          err_msg = 'get_xsqy_tab: failed to read ' // trim(varname) // ' variable'
          call handle_ncerr( nf90_get_var( ncid, varid, xsqy_tab(:,:,:,m) ), trim(err_msg) )
       end do
    endif

    if( .not. is_full_tuv ) then
       del_temp_tab(:ntemp-1) = 1./(temp_tab(2:ntemp) - temp_tab(1:ntemp-1))
       del_conc_tab(:nconc-1) = 1./(conc_tab(2:nconc) - conc_tab(1:nconc-1))
    endif
    dw(:nwave)    = wl(2:nwave+1) - wl(1:nwave)
    w_fac(:nwave) = dw(:nwave)*etfl(:nwave)*1.e-13*wc(:nwave)/hc

    err_msg = 'get_xsqy_tab: failed to close file ' // trim(xsqy_filepath)
    call handle_ncerr( nf90_close( ncid ),trim(err_msg) )

  contains
    subroutine handle_ncerr( ret, mes )
      implicit none
      integer, intent(in) :: ret
      character(len=*), intent(in) :: mes

      if( ret /= nf90_noerr ) then
         write(*,*) 'ERROR: '//trim(mes)//trim(nf90_strerror(ret))
         call exit(ret)
      end if

    end subroutine handle_ncerr


  end subroutine get_xsqy_tab

      subroutine xsqy_int( n, xsqy, tlev, dens_air )
!---------------------------------------------------------------------
!	... interpolate m,t tables for xs * qy
!---------------------------------------------------------------------

      integer, intent(in)  :: n 
      real,    intent(in)  :: tlev(:)
      real,    intent(in)  :: dens_air(:)
      real,    intent(out) :: xsqy(:,:)

      real, parameter :: m0 = 2.45e19
      integer :: tndx, mndx, tndxp1, mndxp1
      integer :: k, ku
      real    :: temp, dens
      real    :: w(4)
      real    :: del_t, del_d

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
        del_t = max( 0.,min( 1.,(temp - temp_tab(tndx))*del_temp_tab(tndx) ) )

!       dens = dens_air(k)
        dens = dens_air(k)/m0
        do mndx = 1,nconc
          if( conc_tab(mndx) > dens ) then
            exit
          endif
        end do
        mndx = max( min( mndx,nconc ) - 1,1 )
        mndxp1 = mndx + 1
        del_d = max( 0.,min( 1.,(dens - conc_tab(mndx))*del_conc_tab(mndx) ) )

        w(1) = (1. - del_t)*(1. - del_d)
        w(2) = del_t*(1. - del_d)
        w(3) = (1. - del_t)*del_d
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

      real,    intent(in)  :: tlev(:)
      real,    intent(in)  :: xs_tab(:,:)
      real,    intent(out) :: xs(:,:)

      integer :: tndx, tndxp1
      integer :: k, ku
      real    :: temp
      real    :: w(2)
      real    :: del_t

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
        del_t = max( 0.,min( 1.,(temp - temp_tab(tndx))*del_temp_tab(tndx) ) )

        w(1) = (1. - del_t)
        w(2) = del_t

        xs(1:nwave,k) = w(1)*xs_tab(1:nwave,tndx) &
                      + w(2)*xs_tab(1:nwave,tndxp1)
      end do

      end subroutine xs_int

end module module_prates_tuv
