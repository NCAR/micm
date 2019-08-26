      MODULE tuv_subs

      use phot_kind_mod, only: dp
      use phot_kind_mod, only: rk => kind_phot

      IMPLICIT none

      private
      public :: tuv_radfld

      CONTAINS

      SUBROUTINE tuv_radfld( nlambda_start, cld_od_opt, cldfrac, nlev, nwave, &
                             zenith, z, albedo, &
                             aircol, o3col, so2col, no2col, &
                             tauaer300, tauaer400, tauaer600, tauaer999, &
                             waer300, waer400, waer600, waer999, &
                             gaer300, gaer400, gaer600, gaer999, &
                             dtaer, omaer, gaer, dtcld, omcld, gcld, &
                             has_aer_ra_feedback, &
                             qll, dobsi, o3_xs, no2_xs, o2_xs, &
                             so2_xs, wmin, wc, tlev, dto2, radfld, efld, &
                             e_dir, e_dn, e_up, &
                             dir_fld, dwn_fld, up_fld, dt_cld, errmsg, errflg )
!-----------------------------------------------------------------------------
!     ... calculate the radiation field
!-----------------------------------------------------------------------------
  
      use molec_ox_xsect,only : molec_ox_xsect_run
      use rad_trans,     only : rtlink
      use phot_util_mod, only : sphers, airmas
      
!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)  :: nlambda_start
      integer, intent(in)  :: nlev
      integer, intent(in)  :: nwave
      integer, intent(in)  :: cld_od_opt
      real(rk), intent(in)  :: zenith
      real(rk), intent(in)  :: dobsi
      real(rk), intent(in)  :: wmin
      real(rk), intent(in)  :: z(:)
      real(rk), intent(in)  :: albedo(:)
      real(rk), intent(in)  :: aircol(:)
      real(rk), intent(in)  :: o3col(:)
      real(rk), intent(in)  :: so2col(:)
      real(rk), intent(in)  :: no2col(:)
      real(rk), intent(in)  :: tauaer300(:)
      real(rk), intent(in)  :: tauaer400(:)
      real(rk), intent(in)  :: tauaer600(:)
      real(rk), intent(in)  :: tauaer999(:)
      real(rk), intent(in)  :: waer300(:)
      real(rk), intent(in)  :: waer400(:)
      real(rk), intent(in)  :: waer600(:)
      real(rk), intent(in)  :: waer999(:)
      real(rk), intent(in)  :: gaer300(:)
      real(rk), intent(in)  :: gaer400(:)
      real(rk), intent(in)  :: gaer600(:)
      real(rk), intent(in)  :: gaer999(:)
      real(rk), intent(in)  :: qll(:)
      real(rk), intent(in)  :: wc(:)
      real(rk), intent(in)  :: tlev(:)
      real(rk), intent(in)  :: cldfrac(:)
      real(rk), intent(in)  :: o2_xs(:)
      real(rk), intent(in)  :: so2_xs(:)
      real(rk), intent(in)  :: o3_xs(:,:)
      real(rk), intent(in)  :: no2_xs(:,:)
      real(rk), intent(in)  :: dto2(:,:)
      real(rk), intent(out) :: radfld(:,:)
      real(rk), intent(out) :: efld(:,:)
      real(rk), intent(inout)  :: dir_fld(:,:), dwn_fld(:,:), up_fld(:,:)
      real(rk), intent(inout)  :: e_dir(:,:), e_dn(:,:), e_up(:,:)
      real(rk), intent(inout)  :: dt_cld(:)
      real(rk), intent(inout)  :: dtaer(:,:), omaer(:,:), gaer(:,:)
      real(rk), intent(inout)  :: dtcld(:,:), omcld(:,:), gcld(:,:)
      logical, intent(in)  :: has_aer_ra_feedback

      character(len=*), intent(out)   :: errmsg
      integer,          intent(out)   :: errflg

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      integer :: wn
      integer :: n_radlev, n_radlevp1
      integer :: nid(0:nlev)
      real(rk) :: dtrl(nlev,nwave)
      real(rk) :: dto3(nlev,nwave)
      real(rk) :: dtso2(nlev,nwave)
      real(rk) :: dtno2(nlev,nwave)
!     real :: dtcld(nlev,nwave)
!     real :: dtaer(nlev,nwave)
      real(rk) :: dtsnw(nlev,nwave)

!     real :: omcld(nlev,nwave)
!     real :: gcld(nlev,nwave)
!     real :: omaer(nlev,nwave)
!     real :: gaer(nlev,nwave)
      real(rk) :: omsnw(nlev,nwave)
      real(rk) :: gsnw(nlev,nwave)

      real(rk) :: edir(nlev+1)
      real(rk) :: edn(nlev+1)
      real(rk) :: eup(nlev+1)
      real(rk) :: fdir(nlev+1)
      real(rk) :: fdn(nlev+1)
      real(rk) :: fup(nlev+1)
      real(rk) :: dsdh(0:nlev,nlev)

      errmsg = ' '
      errflg = 0
      
      n_radlev = size( radfld,dim=2 )
      n_radlevp1 = n_radlev + 1

      do wn = 1,nwave
        omcld(:,wn) = 0._rk
        omaer(:,wn) = 0._rk
        omsnw(:,wn) = 0._rk
        gcld(:,wn)  = 0._rk
        gaer(:,wn)  = 0._rk
        gsnw(:,wn)  = 0._rk
        dtcld(:,wn) = 0._rk
        dtaer(:,wn) = 0._rk
        dtsnw(:,wn) = 0._rk
      end do

      call odrl( wc, aircol, dtrl )
      call odo3( o3col, o3_xs, dto3, dobsi )
      call setso2( so2col, so2_xs, dtso2 )
      call setno2( no2col, no2_xs, dtno2 )
!-------------------------------------------------------------
! aerosol optical depths
!-------------------------------------------------------------
      if( has_aer_ra_feedback ) then
        call setaer( nlambda_start, wc, tauaer300, tauaer400, &
                     tauaer600, tauaer999, waer300, &
                     waer400, waer600, waer999,     &
                     gaer300, gaer400, gaer600,     &
                     gaer999, dtaer, omaer, gaer )
      endif
!-------------------------------------------------------------
! cloud optical depths (cloud water units = g/m3)
!-------------------------------------------------------------
      call setcld( nlambda_start, cld_od_opt, z, qll, cldfrac, &
                   dtcld, omcld, gcld, errmsg, errflg )
      if (errflg .ne. 0) return
      
 !     dt_cld(:n_radlev) = dtcld(2:n_radlevp1,1)

      call sphers( nlev, z, zenith, dsdh, nid )

      do wn = nlambda_start,nwave
        call rtlink( &
           nlev+1, nlev, nwave, &
           wn, albedo(wn), zenith, &
           dsdh, nid, &
           dtrl,  &
           dto3,  &
           dto2, &
           dtso2, &
           dtno2,  &
           dtcld, omcld, gcld, &
           dtaer, omaer, gaer, &
           dtsnw, omsnw, gsnw, &
           edir, edn, eup, fdir, fdn, fup, errmsg, errflg )
        
        radfld(wn,1:n_radlev) = fdir(1:n_radlev) + fdn(1:n_radlev) + fup(1:n_radlev)
        efld(1:n_radlev,wn)    = edir(1:n_radlev) + edn(1:n_radlev) + eup(1:n_radlev)
        dir_fld(1:n_radlev,wn) = fdir(:n_radlev)
        dwn_fld(1:n_radlev,wn) = fdn(1:n_radlev)
        up_fld(1:n_radlev,wn)  = fup(1:n_radlev)
        e_dir(1:n_radlev,wn)   = edir(1:n_radlev)
        e_dn(1:n_radlev,wn)    = edn(1:n_radlev)
        e_up(1:n_radlev,wn)    = eup(1:n_radlev)
        
      end do

      END SUBROUTINE tuv_radfld

      SUBROUTINE odrl( wc, aircol, dtrl )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute Rayleigh optical depths as a function of altitude and wavelength =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  C       - REAL, number of air molecules per cm^2 at each specified    (O)=*
!=            altitude layer                                                 =*
!=  DTRL    - REAL, Rayleigh optical depth at each specified altitude     (O)=*
!=            and each specified wavelength                                  =*
!-----------------------------------------------------------------------------*

!-----------------------------------------------------------------------------*
!     ...dummy arguments
!-----------------------------------------------------------------------------*
      REAL(rk),    intent(in)  :: aircol(:)
      REAL(rk),    intent(in)  :: wc(:)
      REAL(rk),    intent(out) :: dtrl(:,:)

!-----------------------------------------------------------------------------*
!     ...local variables
!-----------------------------------------------------------------------------*
      INTEGER :: nwave, nlyr
      INTEGER :: wn
      REAL(rk)    :: srayl, wmicrn, xx 
      
      nwave = size( wc )
      nlyr  = size( aircol )
!-----------------------------------------------------------------------------*
! compute Rayleigh cross sections and depths:
!-----------------------------------------------------------------------------*
      DO wn = 1,nwave
!-----------------------------------------------------------------------------*
! Rayleigh scattering cross section from WMO 1985 (originally from
! Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
! An empirical formula for its calculation in the homoshpere, Planet.
! Space Sci., 32, 1467-1468, 1984.
!-----------------------------------------------------------------------------*
        wmicrn =  wc(wn)*1.E-3_rk
        IF( wmicrn <= 0.55_rk ) THEN
          xx = 3.6772_rk + 0.389_rk*wmicrn + 0.09426_rk/wmicrn
        ELSE
          xx = 4.04_rk
        ENDIF
        srayl = 4.02e-28_rk/(wmicrn)**xx
!-----------------------------------------------------------------------------*
! alternate (older) expression from
! Frohlich and Shaw, Appl.Opt. v.11, p.1773 (1980).
!-----------------------------------------------------------------------------*
        dtrl(:nlyr,wn) = aircol(:nlyr)*srayl
      END DO

      END SUBROUTINE odrl

      SUBROUTINE odo3( o3col, o3xs, dto3, dobsi )
!-----------------------------------------------------------------------------
!=  NAME:  Optical Depths of O3
!=  PURPOSE:
!=  Compute ozone optical depths as a function of altitude and wavelength
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  O3XS   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (I)
!=           each specified wavelength and altitude
!=  C      - REAL, ozone vertical column increments, molec cm-2, for each (I)
!=           layer
!=  DTO3   - REAL, optical depth due to ozone absorption at each          (O)
!=           specified altitude at each specified wavelength
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      REAL(rk), intent(in)    :: dobsi
      REAL(rk), intent(in)    :: o3col(:)
      REAL(rk), intent(in)    :: o3xs(:,:)
      REAL(rk), intent(inout) :: dto3(:,:)

      INTEGER :: nlyr, nwave
      INTEGER :: wn
      REAL(rk)    :: dob_at_grnd, scale_fac

      nwave = size(o3xs,dim=1)
      nlyr  = size(o3col)

      if( dobsi == 0._rk ) then
!-----------------------------------------------------------------------------
!  no scaling
!-----------------------------------------------------------------------------
        DO wn = 1,nwave
          dto3(:nlyr,wn) = o3col(:nlyr) * o3xs(wn,:nlyr)
        END DO
      else
!-----------------------------------------------------------------------------
!  scale model o3 column to dobsi
!-----------------------------------------------------------------------------
        dob_at_grnd = sum( o3col(:nlyr) )/2.687e16_rk
        scale_fac   = dobsi/dob_at_grnd
        DO wn = 1,nwave
          dto3(:nlyr,wn) = scale_fac * o3col(:nlyr) * o3xs(wn,:nlyr)
        END DO
      endif

      END SUBROUTINE odo3

     SUBROUTINE setso2( colso2, so2_xs, dtso2 )
!-----------------------------------------------------------------------------
!=  PURPOSE:
!=  Set up an altitude profile of SO2 molecules, and corresponding absorption
!=  optical depths.  Subroutine includes a shape-conserving scaling method
!=  that allows scaling of the entire profile to a given overhead SO2
!=  column amount.
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  SO2_XS - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)
!=           each specified wavelength
!=  DTSO2  - REAL, optical depth due to SO2 absorption at each            (O)
!=           specified altitude at each specified wavelength
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      REAL(rk),    intent(in)  :: colso2(:)
      REAL(rk),    intent(in)  :: so2_xs(:)
      REAL(rk),    intent(out) :: dtso2(:,:)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      integer :: nwave, nlyr
      integer :: wn

      nwave = size( so2_xs )
      nlyr  = size( colso2 )

      DO wn = 1,nwave
        dtso2(:nlyr,wn) = colso2(:nlyr)*so2_xs(wn)
      END DO

      END SUBROUTINE setso2

      SUBROUTINE setno2( colno2, no2_xs, dtno2 )
!-----------------------------------------------------------------------------
!=  NAME:  Optical Depths of no2
!=  PURPOSE:
!=  Compute no2 optical depths as a function of altitude and wavelength
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  NO2_XS - REAL, molecular absoprtion cross section (cm^2) of no2 at    (I)
!=           each specified wavelength and altitude
!=  COLNO2 - REAL, no2 vertical column increments, molec cm-2, for each   (I)
!=           layer
!=  DTNO2  - REAL, optical depth due to no2 absorption at each            (O)
!=           specified altitude at each specified wavelength
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      REAL(rk), intent(in)    :: colno2(:)
      REAL(rk), intent(in)    :: no2_xs(:,:)
      REAL(rk), intent(inout) :: dtno2(:,:)

      INTEGER :: nlyr, nwave
      INTEGER :: wn

      nwave = size(no2_xs,dim=1)
      nlyr  = size(colno2)

      DO wn = 1,nwave
        dtno2(:nlyr,wn) = colno2(:nlyr) * no2_xs(wn,:nlyr)
      END DO

      END SUBROUTINE setno2

      subroutine setaer( nlambda_start, wc, tauaer300, tauaer400, &
                         tauaer600, tauaer999,               &
                         waer300, waer400, waer600, waer999, &
                         gaer300, gaer400, gaer600, gaer999, &
                         dtaer, omaer, gaer )
!----------------------------------------------------------------------
! The routine is based on aerosol treatment in module_ra_rrtmg_sw.F
! INPUT: 
! nzlev: number of specified altitude levels in the working grid
! z: specified altitude working grid   
! Aerosol optical properties at 300, 400, 600 and 999 nm. 
!   tauaer300, tauaer400, tauaer600, tauaer999: Layer AODs
!   waer300, waer400, waer600, waer999: Layer SSAs
!   gaer300, gaer400, gaer600, gaer999: Layer asymmetry parameters

! OUTPUT:
! dtaer: Layer AOD at FTUV wavelengths
! omaer: Layer SSA at FTUV wavelengths
! gaer : Layer asymmetry parameters at FTUV wavelengths
!------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: nlambda_start
      real(rk), intent(in)  :: wc(:)
      real(rk), intent(in)  :: tauaer300(:), tauaer400(:),    &
                           tauaer600(:), tauaer999(:)
      real(rk), intent(in)  :: waer300(:), waer400(:),        &
                           waer600(:), waer999(:)
      real(rk), intent(in)  :: gaer300(:), gaer400(:),        &
                           gaer600(:), gaer999(:)
      real(rk), intent(out) :: dtaer(:,:), omaer(:,:), gaer(:,:)

!-----------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------------------------------
      real(rk), parameter :: thresh = 1.e-9_rk
      integer     :: k, wn, nlyr, nwave
      real(rk)        :: ang, slope, wfac

      nlyr =  size(dtaer,dim=1)
      nwave = size(dtaer,dim=2)

wave_loop: &
      do wn = nlambda_start,nwave
        wfac = wc(wn)*1.e-3_rk - .6_rk
        do k = 1,nlyr-1
!-----------------------------------------------------------------------------
! use angstrom exponent to calculate aerosol optical depth; wc is in nm.  
!-----------------------------------------------------------------------------
          if( tauaer300(k) > thresh .and. tauaer999(k) > thresh ) then
            ang = log(tauaer300(k)/tauaer999(k))/log(0.999_rk/0.3_rk)
            dtaer(k,wn) = tauaer400(k)*(0.4_rk/(wc(wn)*1.e-3_rk))**ang
!-----------------------------------------------------------------------------
! ssa - use linear interpolation/extrapolation
!-----------------------------------------------------------------------------
            slope = 5._rk*(waer600(k) - waer400(k))
            omaer(k,wn) = slope*wfac + waer600(k)
            omaer(k,wn) = max( .4_rk,min( 1._rk,omaer(k,wn) ) )
!-----------------------------------------------------------------------------
! asymmetry parameter - use linear interpolation/extrapolation
!-----------------------------------------------------------------------------
            slope = 5._rk*(gaer600(k) - gaer400(k))
            gaer(k,wn) = slope*wfac + gaer600(k)
            gaer(k,wn) = max( .5_rk,min( 1._rk,gaer(k,wn) ) )
          endif
        end do
      end do wave_loop

      end subroutine setaer

      subroutine setcld( nlambda_start, cld_od_opt, z, xlwc, cldfrac, &
                         dtcld, omcld, gcld, errmsg, errflg )
!-----------------------------------------------------------------------------
!= PURPOSE:
!= Set up cloud optical depth, single albedo and g
!-----------------------------------------------------------------------------
!= PARAMETERS:
!= PARAMETERS:
!= NZ - INTEGER, number of specified altitude levels in the working (I)
!= grid
!= Z - real(dp), specified altitude working grid (km) (I)
!= XLWC Cloud water content g/M3 (I)
!=
!= dtcld - cloud optical depth
!= omcld - cloud droplet single albedo
!= gcld  - g
!-----------------------------------------------------------------------------
!
! VERTICAL DOMAIN is from bottom(1) to TOP (TOP=nz)
! CCM from top(1) to bottom(nz)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: nlambda_start
      integer, intent(in) :: cld_od_opt
      real(rk), intent(in)  :: z(:)
      real(rk), intent(in)  :: xlwc(:)
      real(rk), intent(in)  :: cldfrac(:)
      real(rk), intent(inout) :: dtcld(:,:)
      real(rk), intent(inout) :: omcld(:,:)
      real(rk), intent(inout) :: gcld(:,:)

      character(len=*), intent(out)   :: errmsg
      integer,            intent(out)   :: errflg

!-----------------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------------
      real(rk), parameter :: km2m = 1.e3_rk        ! kilometer to meter
      real(rk), parameter :: wden = 1.e6_rk        ! g/m3 (1 m3 water = 1e6 g water)
      real(rk), parameter :: re = 10.0_rk * 1.e-6_rk  ! assuming cloud drop radius = 10 um to M
      real(rk), parameter :: fac = 1._rk/(wden*re)

      integer  :: astat
      integer  :: wn
      integer  :: nlyr, nwave
      real(rk), allocatable :: wrk(:), layer_cldfrac(:)

      errmsg = ' '
      errflg = 0

      nlyr  = size(dtcld,dim=1)
      nwave = size(dtcld,dim=2)

      allocate( wrk(nlyr),layer_cldfrac(nlyr),stat=astat )
      if( astat /= 0 ) then
         errmsg = 'setcld: failed to allocate wrk'
         errflg = 1
         return
      endif

!-----------------------------------------------------------------------------
! ... calculate optical depth
!-----------------------------------------------------------------------------     
      wrk(1:nlyr-1) = (z(2:nlyr) - z(1:nlyr-1))*km2m   !  (km -> m)
      wrk(1:nlyr-1) = 1.5_rk * .5_rk*(xlwc(1:nlyr-1) + xlwc(2:nlyr))*wrk(1:nlyr-1)*fac
      wrk(1:nlyr-1) = max( wrk(1:nlyr-1),0._rk )
      if( cld_od_opt == 2 ) then
        layer_cldfrac(1:nlyr-1) = .5_rk*(cldfrac(1:nlyr-1) + cldfrac(2:nlyr))
        wrk(1:nlyr-1) = wrk(1:nlyr-1)*layer_cldfrac(1:nlyr-1)*sqrt( layer_cldfrac(1:nlyr-1) )
      endif
!----------------------------------------------------
! ....calculate cloud optical depth T
! following Liao et al. JGR, 104, 23697, 1999
!----------------------------------------------------
      if( any( wrk(1:nlyr-1) > 0._rk ) ) then
        do wn = nlambda_start,nwave
          dtcld(1:nlyr-1,wn) = wrk(1:nlyr-1)
          omcld(1:nlyr-1,wn) = .9999_rk
          gcld (1:nlyr-1,wn) = .85_rk
        end do
      endif

      if( allocated( wrk ) ) then
        deallocate( wrk )
      endif
      if( allocated( layer_cldfrac ) ) then
        deallocate( layer_cldfrac )
      endif

      end subroutine setcld


      END MODULE tuv_subs
