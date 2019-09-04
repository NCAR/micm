module  module_prates_tuv
  use phot_kind_mod, only: rk => kind_phot
  use module_rxn, only : xsqy_table => xsqy_tab, the_subs, npht_tab, rxn_init
  use module_rxn, only : get_initialization, set_initialization
  use wavelength_grid,only: nwave, wl, wc
  use params_mod, only : hc, qnan
  
  implicit none

  private
  public :: calc_tuv_init
  public :: calc_tuv_prates
  public :: rxn_ndx
  public :: read_etf
  public :: update_etf

  integer :: nconc, ntemp

  real(rk), allocatable :: xsqy_tab(:,:,:,:)
  real(rk), allocatable :: dw(:)
  real(rk), allocatable :: photon_flux(:) ! /cm2/sec
  real(rk), allocatable :: temp_tab(:)
  real(rk), allocatable :: conc_tab(:)
  real(rk), allocatable :: del_temp_tab(:)
  real(rk), allocatable :: del_conc_tab(:)
  character(len=32), allocatable :: tuv_jname(:)

  integer :: nj 
  integer :: j_o2_ndx = -1
  integer :: jo2_a_ndx = -1
  integer :: jo2_b_ndx = -1
  logical :: is_full_tuv = .true.

  logical, allocatable :: xsqy_is_zdep(:)
  integer, protected, allocatable :: rxn_ndx(:)

  real(rk) :: esfact = 1.0_rk
  real(rk) :: xnan

contains

  !------------------------------------------------------------------------------
  ! initialize
  !------------------------------------------------------------------------------
  subroutine calc_tuv_init( full_tuv, nlevs, jnames, xsqy_filepath, errmsg, errflg )

    logical, intent(in) :: full_tuv
    integer, intent(in) :: nlevs
    character(len=*), intent(in) :: jnames(:)
    character(len=*), intent(in) :: xsqy_filepath
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    
    integer :: i,j,n,jndx
    integer :: astat
    logical :: rxn_initialized
    real(rk) :: dummy(nlevs)

    xnan = qnan()
    
    is_full_tuv = full_tuv
    nj = size(jnames)
    allocate( tuv_jname(nj), stat=astat )
    if( astat /= 0 ) then
       errmsg = 'calc_tuv_init: failed to allocate tuv_jname'
       errflg = astat
       return
    endif

    tuv_jname(:) = jnames(:)

    find_jo2: do n=1,nj
       if (trim(tuv_jname(n)) == 'j_o2') j_o2_ndx = n
       if (trim(tuv_jname(n)) == 'jo2_a') jo2_a_ndx = n
       if (trim(tuv_jname(n)) == 'jo2_b') jo2_b_ndx = n
    end do find_jo2

    if( is_full_tuv ) then

       call rxn_init( nwave+1, wl, errmsg, errflg )
       allocate( rxn_ndx(nj), stat=astat )
       if( astat /= 0 ) then
          errmsg = 'calc_tuv_init: failed to allocate rxn_ndx'
          errflg = astat
          return
       endif
       rxn_ndx(1:nj) = -1
       do j = 1,nj
          if( j == j_o2_ndx ) then
             rxn_ndx(j) = 0
          else
             do n = 2,npht_tab
                if( trim(xsqy_table(n)%rxn_name) == trim(tuv_jname(j)) ) then
                   if (.not.any(rxn_ndx==n)) then
                      rxn_ndx(j) = n
                      exit
                   else
                      errmsg = trim(errmsg)//' '//trim(tuv_jname(j))
                      errflg = 1
                   end if
                endif
             enddo
          endif
       enddo

       if (errflg/=0) then
          errmsg = 'calc_tuv_init: deplicate jnames: '//trim(errmsg)
          return
       end if
       
       if (any(rxn_ndx(1:nj)<0)) then
          errflg = 1
          errmsg = 'calc_tuv_init: not recognized jnames:'
          do j = 1,nj
             if (rxn_ndx(j) < 0 ) then
                errmsg = trim(errmsg)//' '//trim(tuv_jname(j))
             end if
          end do
          return
       end if

       rxn_initialized = .not. get_initialization()
       if( .not. rxn_initialized ) then
          do n = 1,nj
             jndx = rxn_ndx(n)
             if( jndx > 0 ) then
                call the_subs(jndx)%xsqy_sub(nwave+1,wl,wc,nlevs,dummy,dummy,jndx, errmsg, errflg )
             endif
          enddo
          call set_initialization( status=.false. )
       endif

    else
       call read_xsqy_tables(xsqy_filepath, errmsg, errflg)
       allocate( xsqy_is_zdep(nj), stat=astat )
       if( astat /= 0 ) then
          errmsg = 'calc_tuv_init: failed to allocate xsqy_is_zdep'
          errflg = astat
          return
       endif
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
    endif

  end subroutine calc_tuv_init

  !------------------------------------------------------------------------------
  ! compute Js for a column given photon fluxes in each grid box of the column (nlevs)
  !------------------------------------------------------------------------------
  subroutine calc_tuv_prates(kts,kte,nlevs, tlev, dens_air, rad_fld, srb_o2_xs, tuv_prate, errmsg, errflg)

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
    real(rk) :: sq2d(nlevs,nwave)
    real(rk) :: sq1d(nwave,1)

    tuv_prate = xnan
    rad_fld_tpose = xnan

    if( .not. is_full_tuv ) then
       if( any( .not. xsqy_is_zdep(:) ) ) then
          rad_fld_tpose = transpose( rad_fld )
       endif
    elseif( any( xsqy_table(1:nj)%tpflag == 0 ) ) then
       rad_fld_tpose = transpose( rad_fld )
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
                   call the_subs(jndx)%xsqy_sub(nwave+1,wl,wc,nlevs,tlev,dens_air,jndx, errmsg, errflg, sq=sq2d)
                else
                   call the_subs(jndx)%xsqy_sub(nwave+1,wl,wc,nlevs,tlev,dens_air,jndx, errmsg, errflg, sq=sq1d)
                end if
             endif
          endif
       else
          if( is_full_tuv ) then
             sq1d = 1._rk
          else
             call sjo2( kte, nwave, srb_o2_xs, xsqy )
          endif
       endif
       !---------------------------------------------------------------------
       ! compute tuv photorates
       !---------------------------------------------------------------------
       if( .not. is_full_tuv ) then
          if( xsqy_is_zdep(n) ) then
             do k = kts,kte
                xsect(1:nwave) = xsqy(1:nwave,k)*photon_flux(1:nwave)*esfact
                tuv_prate(k,n) = dot_product( rad_fld(1:nwave,k),xsect(1:nwave) )
             end do
          else
             xsect(1:nwave) = xsqy_tab(1:nwave,1,1,n)*photon_flux(1:nwave)*esfact
             tuv_prate(:,n) = matmul( rad_fld_tpose(:,1:nwave),xsect(1:nwave) )
          endif
       else
          if(n==jo2_a_ndx .or. n==jo2_b_ndx .or. n==j_o2_ndx) then
             do k = kts,kte
                xsect(1:nwave) = sq1d(1:nwave,1)*srb_o2_xs(1:nwave,k)*photon_flux(1:nwave)*esfact
                tuv_prate(k,n) = dot_product( rad_fld(1:nwave,k),xsect(1:nwave) )
             end do
          else if ( xsqy_table(jndx)%tpflag > 0 ) then
             do k = kts,kte
                xsect(1:nwave) = sq2d(k,1:nwave)*photon_flux(1:nwave)*esfact
                tuv_prate(k,n) = dot_product( rad_fld(1:nwave,k),xsect(1:nwave) )
             end do
          else
             xsect(1:nwave) = sq1d(1:nwave,1)*photon_flux(1:nwave)*esfact
             tuv_prate(:,n) = matmul( rad_fld_tpose(:,1:nwave),xsect(1:nwave) )
          end if
       endif
    end do rate_loop

  end subroutine calc_tuv_prates

  ! update extraterrestrial photon flux
  subroutine update_etf( etf_in )
    real(rk), intent(in) :: etf_in(nwave) ! watts/m2/nm

    ! nm * (watts/m2/nm) * nw * 1.e-13/(watts m sec) --> /cm2/sec
    photon_flux(:nwave) = dw(:nwave)*etf_in(:nwave)*wc(:nwave)*1.e-13_rk/hc

  end subroutine update_etf
  
  subroutine read_etf(filepath, errmsg, errflg )
    !---------------------------------------------------------------------
    !	... read in ETF
    !---------------------------------------------------------------------
    use netcdf

    !---------------------------------------------------------------------
    !	... dummy args
    !---------------------------------------------------------------------
    character(len=*), intent(in)  :: filepath
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !---------------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------------
    integer :: astat, ret
    integer :: m
    integer :: ncid, dimid, varid
    character(len=64) :: varname

    real(rk) :: etfl(nwave)

    !---------------------------------------------------------------------
    !	... open file
    !---------------------------------------------------------------------
    ret = nf90_open( trim(filepath), nf90_noclobber, ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_etf: failed to open file ' // trim(filepath)
       return
    end if

    !---------------------------------------------------------------------
    !	... allocate module arrays
    !---------------------------------------------------------------------
    allocate( dw(nwave), photon_flux(nwave), stat=astat )
    if( astat /= 0 ) then
       errmsg = 'read_etf: failed to allocate'
       errflg = astat
       return
    endif

    dw(:nwave) = wl(2:nwave+1) - wl(1:nwave)
    photon_flux(:nwave) = xnan

    !---------------------------------------------------------------------
    !	... read arrays
    !---------------------------------------------------------------------
    ret = nf90_inq_varid( ncid, 'etf', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_etf: failed to get etfl variable id'
       return
    end if

    if (varid>0) then
       ret = nf90_get_var( ncid, varid, etfl )
       if( ret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'read_etf: failed to read etfl variable'
          return
       end if

       photon_flux(:nwave) = dw(:nwave)*etfl(:nwave)*1.e-13_rk*wc(:nwave)/hc !  watts/m2/nm --> /cm2/sec
    end if

    !---------------------------------------------------------------------
    !	... close file
    !---------------------------------------------------------------------
    ret = nf90_close( ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_etf: failed to close file ' // trim(filepath)
       return
    end if

  end subroutine read_etf

  subroutine read_xsqy_tables( xsqy_filepath, errmsg, errflg )
    use netcdf

    character(len=*), intent(in)  :: xsqy_filepath
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !---------------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------------
    integer :: astat, ret
    integer :: m
    integer :: ncid, dimid, varid
    character(len=64) :: varname

    errmsg = ' '
    errflg = 0

    ret = nf90_open( trim(xsqy_filepath), nf90_noclobber, ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to open file ' // trim(xsqy_filepath)
       return
    end if

    ret = nf90_inq_dimid( ncid, 'ntemp', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to get ntemp id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=ntemp )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to get ntemp'
       return
    end if
    ret = nf90_inq_dimid( ncid, 'nconc', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to get nconc id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=nconc )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to get nconc'
       return
    end if

    allocate( temp_tab(ntemp), conc_tab(nconc), stat=astat )
    if( astat /= 0 ) then
       errmsg = 'read_xsqy_tables: failed to allocate'
       errflg = astat
       return
    endif
    allocate( del_temp_tab(ntemp-1), del_conc_tab(nconc-1), stat=astat )
    if( astat /= 0 ) then
       errmsg = 'read_xsqy_tables: failed to allocate'
       errflg = astat
       return
    endif
    allocate( xsqy_tab(nwave,ntemp,nconc,nj), stat=astat )
    if( astat /= 0 ) then
       errmsg = 'read_xsqy_tables: failed to allocate'
       errflg = astat
       return
    endif

    ret = nf90_inq_varid( ncid, 'temps', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to temp_tab variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, temp_tab )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to read temp_tab variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'concs', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to conc_tab variable id'
       return
    end if

    ret = nf90_get_var( ncid, varid, conc_tab )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to read conc_tab variable'
       return
    end if
    do m = 1,nj
       varname = trim(tuv_jname(m)) // '_xsqy'
       ret = nf90_inq_varid( ncid, trim(varname), varid )
       if( ret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'read_xsqy_tables: failed to ' // trim(varname) //' variable id'
          return
       end if
       ret = nf90_get_var( ncid, varid, xsqy_tab(:,:,:,m) )
       if( ret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'read_xsqy_tables: failed to read ' // trim(varname) // ' variable'
          return
       end if
    end do

    del_temp_tab(:ntemp-1) = 1._rk/(temp_tab(2:ntemp) - temp_tab(1:ntemp-1))
    del_conc_tab(:nconc-1) = 1._rk/(conc_tab(2:nconc) - conc_tab(1:nconc-1))

    ret = nf90_close( ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_xsqy_tables: failed to close file ' // trim(xsqy_filepath)
       return
    end if
    
  end subroutine read_xsqy_tables

 
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
