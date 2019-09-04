!=============================================================================
! This file contains the following subroutines, related to the calculation
! of radiation at Lyman-alpha and Schumann-Runge wavelengths:
!     la_srb
!     lymana
!     schum
!     effxs
!     calc_params
!     init_xs
! and the following functions
!     chebev
!=============================================================================

      module la_srb_mod

      use phot_kind_mod, only: DP
      use phot_kind_mod, only: rk => kind_phot

      implicit none

      private

      public :: la_srb_comp
      public :: la_srb_init

      integer, parameter :: kla = 2
      integer, parameter :: ksrb = 18
      integer, parameter :: nla =  kla - 1
      integer, parameter :: nsrb = ksrb - 1

      integer :: nchebev_term=-1, nchebev_wave=-1

      integer :: ila, isrb
      real(kind=dp), allocatable :: chebev_ac(:,:)
      real(kind=dp), allocatable :: chebev_bc(:,:)

      ! Lyman-Alpha wavelength band edges
      real(rk), parameter :: wlla(kla) = (/ 121.4_rk, 121.9_rk /)

      ! Schumann-Runge wavelength band edges
      real(rk), parameter :: wlsrb(ksrb) = &
           (/ 174.4_rk, 177.0_rk, 178.6_rk, 180.2_rk, 181.8_rk, &
              183.5_rk, 185.2_rk, 186.9_rk, 188.7_rk, 190.5_rk, &
              192.3_rk, 194.2_rk, 196.1_rk, 198.0_rk, 200.0_rk, &
              202.0_rk, 204.1_rk, 205.8_rk/) ! 17 SRB bands

      real(rk) :: xnan

      contains

      subroutine la_srb_init( errmsg, errflg )
        use params_mod, only: input_data_root, qnan
        use wavelength_grid, only: nwave, wc
        use netcdf

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        integer :: ncid, dimid, varid, iw, i
        integer :: astat, ret
        character(len=512) :: filepath

        xnan = qnan()

        filepath = trim(input_data_root)//'/chebev_coeffs.nc'

        errmsg = ' '
        errflg = 0

        ! open file
        ret = nf90_open( trim(filepath), nf90_noclobber, ncid )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to open '//trim(filepath)
           return
        end if

        ret = nf90_inq_dimid( ncid, 'nchebev_term', dimid )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to get nchebev_term id'
           return
        end if
        ret = nf90_inquire_dimension( ncid, dimid, len=nchebev_term )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to get nchebev'
           return
        end if
        ret = nf90_inq_dimid( ncid, 'nchebev_wave', dimid )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to get nchebev_wave id'
           return
        end if
        ret = nf90_inquire_dimension( ncid, dimid, len=nchebev_wave )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to get nchebev'
           return
        end if

        allocate( chebev_ac(nchebev_term,nchebev_wave), chebev_bc(nchebev_term,nchebev_wave), stat=astat )

        if( astat /= 0 ) then
           errflg = astat
           errmsg = 'la_srb_init: failed to allocate chebev memory'
           return
        end if
        ret = nf90_inq_varid( ncid, 'chebev_ac', varid )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to get chebev_ac variable id'
           return
        end if
        ret = nf90_get_var( ncid, varid, chebev_ac )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to read chebev_ac variable'
           return
        end if
        ret = nf90_inq_varid( ncid, 'chebev_bc', varid )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to get chebev_bc variable id'
           return
        end if
        ret = nf90_get_var( ncid, varid, chebev_bc )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to read chebev_bc variable'
           return
        end if

        ! close the file
        ret = nf90_close( ncid )
        if( ret /= nf90_noerr ) then
           errflg = 1
           errmsg = 'la_srb_init: failed to close '//trim(filepath)
           return
        end if

        ! check that the wavelength grid includes Lyman-alpha and Schumann-Runge wavelength bands

        ila = -1
        isrb = -1
     
        ila_loop: do iw = 1,nwave-1
           if (wc(iw)>wlla(1) .and. wc(iw)<wlla(2) .and. wc(iw+1)>wlla(2)) then
              ila = iw
              exit ila_loop
           end if
        end do ila_loop

       isrb_loop: do iw = 1,nwave-nsrb
           if (wc(iw)>wlsrb(1) .and. wc(iw)<wlsrb(2)) then
              do i = 1,nsrb-1
                 if ( .not. (wc(iw+i)>wlsrb(i+1) .and. wc(iw+i)<wlsrb(i+2)) ) then
                    exit isrb_loop
                 endif
              end do
              if ( .not. (wc(iw+nsrb)>wlsrb(nsrb+1)) ) then
                 exit isrb_loop
              end if                                  
              isrb = iw
              exit isrb_loop
           end if
        end do isrb_loop

        if (ila<1 .or. isrb<1) then
           errflg = 1
           errmsg = 'la_srb_init: wavelength grid must contain Lyman-alpha and Schumann-Runge wavelength bands'
           return
        end if

      END SUBROUTINE la_srb_init

      SUBROUTINE la_srb_comp( nlyr, wmin, tlev, vcol, scol, o2vmr, o2_xs, dto2, srb_o2_xs, errmsg, errflg )
!-----------------------------------------------------------------------------
!=  PURPOSE:
!=  Compute equivalent optical depths for O2 absorption, and O2 effective
!=  absorption cross sections, parameterized in the Lyman-alpha and SR bands
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)
!=            grid
!=  Z       - REAL, specified altitude working grid (km)                  (I)
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)
!=            wavelength grid
!=  WL      - REAL, vector of lxower limits of wavelength intervals in    (I)
!=            working wavelength grid
!=  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)
!=            altitude layer
!=  ZEN     - REAL, solar zenith angle                                    (I)
!=
!=  O2XS1   - REAL, O2 cross section from rdo2xs                          (I)
!=
!=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)
!=            vertical layer at each specified wavelength
!=  O2XS    - REAL, molecular absorption cross section in SR bands at     (O)
!=            each specified altitude and wavelength.  Includes Herzberg
!=            continuum.
!-----------------------------------------------------------------------------

      use params_mod, only : largest

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      INTEGER, intent(in) :: nlyr
      REAL(rk), intent(in) :: wmin
      REAL(rk), intent(in) :: tlev(:)

      REAL(rk), intent(in) :: vcol(:)
      REAL(rk), intent(in) :: scol(:)
      REAL(rk), intent(in) :: o2vmr(:)
      REAL(rk), intent(in) :: o2_xs(:)
      REAL(rk), intent(inout) :: dto2(:,:)
      REAL(rk), intent(inout) :: srb_o2_xs(:,:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      REAL(rk) :: secchi(nlyr)
      REAL(rk) :: o2col(nlyr)

!-----------------------------------------------------------------------------
! Lyman-alpha variables
! O2 optical depth and equivalent cross section in the Lyman-alpha region
!-----------------------------------------------------------------------------
!     INTEGER :: nlev
      INTEGER :: nlev_srb
      INTEGER :: k, iw, wn
      REAL(rk)    :: dto2la(nlyr,nla), o2xsla(nlyr,nla)

!-----------------------------------------------------------------------------
! grid on which Koppers' parameterization is defined
! O2 optical depth and equivalent cross section on Koppers' grid
!-----------------------------------------------------------------------------
      REAL(rk)    :: dto2k(nlyr,nsrb), o2xsk(nlyr,nsrb)

      errmsg = ' '
      errflg = 0

     if (nchebev_wave<1) then
         errflg = 1
         errmsg = 'la_srb_comp not initialized'
         srb_o2_xs = xnan
         dto2 = xnan
         return
      end if
      
      nlev_srb = size( srb_o2_xs,dim=2 )
!----------------------------------------------------------------------
! initalize O2 cross sections 
!----------------------------------------------------------------------
      DO k = 1, nlev_srb
        srb_o2_xs(:,k) = o2_xs(:)
      END DO

      IF( wmin <= wlsrb(nsrb) ) THEN
!----------------------------------------------------------------------
! Slant O2 column and x-sections.
!----------------------------------------------------------------------
        o2col(:nlyr) = o2vmr(:nlyr) * scol(:nlyr)
!----------------------------------------------------------------------
! Effective secant of solar zenith angle.  
! Use 2.0 if no direct sun (value for isotropic radiation)
! For nz, use value at nz-1
!----------------------------------------------------------------------
        WHERE( scol(:nlyr) > .1_rk*largest ) 
          secchi(:nlyr) = 2._rk
        ELSEWHERE
          secchi(:nlyr) = scol(:nlyr)/vcol(:nlyr)
        ENDWHERE

!---------------------------------------------------------------------
! Lyman-Alpha parameterization, output values of O2 optical depth
! and O2 effective (equivalent) cross section
!----------------------------------------------------------------------
        CALL lymana( nlyr, o2col, secchi, dto2la, o2xsla )
        DO wn = ila, ila + nla - 1
          iw = wn - ila + 1
          dto2(:nlyr,wn)          = dto2la(:nlyr,iw) 
          srb_o2_xs(wn, 1:nlyr) = o2xsla(1:nlyr,iw)
          srb_o2_xs(wn,nlyr+1:nlev_srb) = o2xsla(nlyr,iw)
        ENDDO

!------------------------------------------------------------------------------
! Koppers' parameterization of the SR bands, output values of O2
! optical depth and O2 equivalent cross section 
!------------------------------------------------------------------------------
        CALL schum( nlyr, o2col, tlev, secchi, dto2k, o2xsk )
        DO wn = isrb, isrb + nsrb - 1
          iw = wn - isrb + 1
          dto2(:nlyr,wn)          = dto2k(:nlyr,iw)
          srb_o2_xs(wn, 1:nlyr)   = o2xsk(1:nlyr,iw)
          srb_o2_xs(wn,nlyr+1:nlev_srb) = o2xsk(nlyr,iw)
        ENDDO
      ENDIF

      END SUBROUTINE la_srb_comp

      SUBROUTINE lymana( nlyr, o2col, secchi, dto2la, o2xsla )
!-----------------------------------------------------------------------------
!=  PURPOSE:
!=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha
!=  bands and an effective O2 optical depth at all altitudes.  Parameterized
!=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the
!=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,
!=  Vol.24, No.21, pp 2659-2662, 1997.
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)
!=            grid
!=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)
!=            altitude
!=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)
!=            vertical layer
!=  O2XSLA  - REAL, molecular absorption cross section in LA bands        (O)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      INTEGER, intent(in) :: nlyr
      REAL(rk),    intent(in) :: o2col(:)
      REAL(rk),    intent(in) :: secchi(:)
      REAL(rk), intent(inout) :: dto2la(nlyr,nla), o2xsla(nlyr,nla)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      REAL(rk), parameter    :: xsmin = 1.e-20_rk
      REAL(kind=DP), parameter :: rmmin = 1.e-100_DP

      INTEGER :: k, kp1, wn
      REAL(kind=DP) :: o2_col
      REAL(kind=DP) :: rm(nlyr), ro2(nlyr)
      REAL(kind=DP) :: rm_wrk(3), ro2_wrk(3)

      real(kind=dp), parameter :: b(3) = (/ 6.8431e-01_DP,  2.29841e-01_DP,  8.65412e-02_DP /)
      real(kind=dp), parameter :: c(3) = (/ 8.22114e-21_DP, 1.77556e-20_DP,  8.22112e-21_DP /)
      real(kind=dp), parameter :: d(3) = (/ 6.0073e-21_DP,  4.28569e-21_DP,  1.28059e-20_DP /)
      real(kind=dp), parameter :: e(3) = (/ 8.21666e-21_DP, 1.63296e-20_DP,  4.85121e-17_DP /)

      do wn = 1,nla
        dto2la(:nlyr,wn) = 0._rk
        o2xsla(:nlyr,wn) = 0._rk
      end do
!-----------------------------------------------------------------------------
! calculate reduction factors at every layer
!-----------------------------------------------------------------------------
      rm(:nlyr)  = 0._DP
      ro2(:nlyr) = 0._DP
      DO k = 1, nlyr
        o2_col = real( o2col(k),DP )
        rm_wrk(:)  = b(:) * EXP( -c(:) * o2_col )
        ro2_wrk(:) = d(:) * EXP( -e(:) * o2_col )
        rm(k)  = sum( rm_wrk )
        ro2(k) = sum( ro2_wrk )
      ENDDO

!-----------------------------------------------------------------------------
! calculate effective O2 optical depths and effective O2 cross sections
!-----------------------------------------------------------------------------
      DO k = 1, nlyr-1
        kp1 = k + 1
        IF (rm(k) > rmmin) THEN
          IF (ro2(k) > rmmin) THEN
            o2xsla(k,1) = REAL( ro2(k)/rm(k) ,rk)
          ELSE
            o2xsla(k,1) = xsmin
          ENDIF

          IF (rm(kp1) > 0._DP) THEN
            dto2la(k,1) = LOG( rm(kp1) )/secchi(kp1)  &
                        - LOG( rm(k))   /secchi(k)
          ELSE
            dto2la(k,1) = 1000._rk
          ENDIF
        ELSE
          dto2la(k,1) = 1000._rk
          o2xsla(k,1) = xsmin
        ENDIF
      END DO

!-----------------------------------------------------------------------------
! do top layer separately
!-----------------------------------------------------------------------------
      IF( rm(nlyr) > rmmin ) THEN
        o2xsla(nlyr,1) = REAL( ro2(nlyr)/rm(nlyr) ,rk)
      ELSE
        o2xsla(nlyr,1) = xsmin
      ENDIF

      END SUBROUTINE lymana

      SUBROUTINE schum( nlyr, o2col, tlev, secchi, dto2, o2xsk )
!-----------------------------------------------------------------------------
!=  PURPOSE:
!=  Calculate the equivalent absorption cross section of O2 in the SR bands.
!=  The algorithm is based on parameterization of G.A. Koppers, and
!=  D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]
!=  Final values do include effects from the Herzberg continuum.
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)
!=            grid
!=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)
!=            altitude
!=  TLEV    - tmeperature at each level
!=  SECCHI  - ratio of slant to vertical o2 columns
!=  DTO2    - REAL, optical depth due to O2 absorption at each specified
!=            vertical layer at each specified wavelength
!=  O2XSK  - REAL, molecular absorption cross section in SR bands at
!=            each specified wavelength.  Includes Herzberg continuum
!-----------------------------------------------------------------------------

      use params_mod, only : precis

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      INTEGER, intent(in) :: nlyr
      REAL(rk),    intent(in) :: o2col(:)
      REAL(rk),    intent(in) :: tlev(:), secchi(:)
      REAL(rk), intent(inout) :: dto2(:,:), o2xsk(:,:)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      REAL(rk), parameter :: o2col_min = exp( 38._rk )

      INTEGER :: wn, k, ktop, ktop1, kbot, nlyrm1
      REAL(rk)    :: x
      REAL(rk)    :: o2col1(nlyr)
      REAL(rk)    :: xs(nsrb)
      real(rk), parameter :: xslod(nsrb) = &
           (/ 6.2180730E-21_rk, 5.8473627E-22_rk, 5.6996334E-22_rk, &
              4.5627094E-22_rk, 1.7668250E-22_rk, 1.1178808E-22_rk, &
              1.2040544E-22_rk, 4.0994668E-23_rk, 1.8450616E-23_rk, &
              1.5639540E-23_rk, 8.7961075E-24_rk, 7.6475608E-24_rk, &
              7.6260556E-24_rk, 7.5565696E-24_rk, 7.6334338E-24_rk, &
              7.4371992E-24_rk, 7.3642966E-24_rk /)

      nlyrm1 = nlyr - 1
!-----------------------------------------------------------------------------
!     ...Initialize cross sections to values at large optical depth
!-----------------------------------------------------------------------------
      DO wn = 1, nsrb
        o2xsk(:nlyr,wn) = xslod(wn)
      END DO

!-----------------------------------------------------------------------------
!     Calculate cross sections
!     Set smallest O2col = exp(38.) molec cm-2
!     to stay in range of parameterization
!     given by Koppers et al. at top of atm.
!-----------------------------------------------------------------------------
      ktop = 2*nlyr
      kbot = 0

      DO k = 1,nlyr
        o2col1(k) = MAX( o2col(k),o2col_min )
        x  = LOG( o2col1(k) )
        IF (x < 38.0_rk) THEN
          ktop1 = k-1
          ktop  = MIN(ktop1,ktop)
        ELSE IF (x > 56.0_rk) THEN
          kbot = k
        ELSE
          CALL effxs( x, tlev(k), xs )
          o2xsk(k,:nsrb) = xs(:nsrb)
        ENDIF
      END DO

!-----------------------------------------------------------------------------
!  fill in cross section where X is out of range by repeating edge table values
!  Do not allow kbot = nlyr to avoid division by zero in no light case.
!-----------------------------------------------------------------------------
      IF( kbot == nlyr) then
        kbot = nlyrm1
      ENDIF

      IF( kbot > 0 ) THEN
        DO wn = 1,nsrb
          o2xsk(:kbot,wn) = o2xsk(kbot+1,wn)
        END DO
      ENDIF

      IF( ktop < nlyr ) THEN
        DO wn = 1,nsrb
          o2xsk(ktop+1:nlyr,wn) = o2xsk(ktop,wn)
        END DO
      ENDIF

!-----------------------------------------------------------------------------
!  Calculate incremental optical depths 
!-----------------------------------------------------------------------------
      dto2(nlyr,1:nsrb) = 0.0_rk       ! set optical depth to zero at top
      DO wn = 1,nsrb
!-----------------------------------------------------------------------------
!     ... calculate an optical depth weighted by density,
!         put in mean value estimate, if in shade
!-----------------------------------------------------------------------------
        WHERE (ABS(1._rk - o2col1(2:nlyr)/o2col1(:nlyrm1)) <= 2._rk*precis)
          dto2(:nlyrm1,wn) = o2xsk(2:nlyr,wn)*o2col1(2:nlyr)/real(nlyrm1)
        ELSEWHERE
          dto2(:nlyr-1,wn) = ABS( &
            (o2xsk(2:nlyr,wn)*o2col1(2:nlyr) - o2xsk(:nlyrm1,wn)*o2col1(:nlyrm1)) &
            /(1._rk + LOG(o2xsk(2:nlyr,wn)/o2xsk(:nlyrm1,wn))  &
              / LOG(o2col1(2:nlyr)/o2col1(:nlyrm1))) )
!-----------------------------------------------------------------------------
!     ... change to vertical optical depth
!-----------------------------------------------------------------------------
          dto2(:nlyrm1,wn) = 2._rk * dto2(:nlyrm1,wn)/(secchi(:nlyr-1)+secchi(2:nlyr))
        ENDWHERE
      END DO 

      END SUBROUTINE schum

      SUBROUTINE EFFXS( x, t, xs )
!-----------------------------------------------------------------------------
!     Subroutine for evaluating the effective cross section
!     of O2 in the Schumann-Runge bands using parameterization
!     of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
!     68-79, 1996]
!      
!     method:
!     ln(xs) = A(X)[T-220]+B(X)
!     X = log of slant column of O2
!     A,B calculated from Chebyshev polynomial coeffs
!     AC and BC using NR routine chebev.  Assume interval
!     is 38<ln(NO2)<56.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      REAL(rk), intent(in)  :: t, x
      REAL(rk), intent(out) :: xs(nsrb)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      REAL(rk)    :: a(nsrb), b(nsrb) 

      call calc_params( x, a, b )

      xs(:nsrb) = EXP( a(:nsrb)*( t - 220._rk) + b(:nsrb) )

      END SUBROUTINE EFFXS

      SUBROUTINE CALC_PARAMS( x, a, b )
!-----------------------------------------------------------------------------
!     calculates coefficients (A,B), used in calculating the
!     effective cross section, for nsrb wavelength intervals
!     as a function of log O2 column density (X)
!     Wavelength intervals are defined in WMO1985
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      REAL(rk), intent(in)  :: x
      REAL(rk), intent(out) :: a(nsrb), b(nsrb)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      INTEGER :: wn

!-----------------------------------------------------------------------------
!     call Chebyshev Evaluation routine to calc A and B from
!     set of 20 coeficients for each wavelength
!-----------------------------------------------------------------------------

      DO wn = 1,nsrb
        a(wn) = chebev( 38.0_rk , 56.0_rk, chebev_ac(:,wn), nchebev_term, x )
        b(wn) = chebev( 38.0_rk , 56.0_rk, chebev_bc(:,wn), nchebev_term, x )
      END DO

      END SUBROUTINE CALC_PARAMS

      REAL(rk) FUNCTION chebev( a, b, c, m, x )
!-------------------------------------------------------------
!     Chebyshev evaluation algorithm
!     See Numerical recipes p193
!-------------------------------------------------------------
      
!-------------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------------
      INTEGER, intent(in) :: m
      REAL(rk),    intent(in) :: a, b, x
      REAL(kind=DP), intent(in) :: c(:)

!-------------------------------------------------------------
!       ... local variables
!-------------------------------------------------------------
      INTEGER :: j
      REAL(rk)    :: d, dd, sv, y, y2

      IF( (x - a)*(x - b) > 0._rk) THEN
	chebev = 0.0_rk
      ELSE
	d  = 0._rk
        dd = 0._rk
        y  = (2._rk*x - a - b)/(b - a)
        y2 = 2._rk*y
        DO J = m,2,-1
          sv = d
          d  = y2*d - dd + real( c(J) )
          dd = sv
        END DO
        chebev = y*d - dd + 0.5_rk*real( c(1) )
      ENDIF
	
      END FUNCTION chebev

      end module la_srb_mod
