! This file contains subroutines related to reading the
! absorption cross sections of gases that contribute to atmospheric transmission:
! Some of these subroutines are also called from rxn.f when loading photolysis cross sections
! for these same gases. It is possible to have different cross sections for 
! transmission and for photolysis, e.g. for ozone, Bass et al. could be used
! for transmission while Molina and Molina could be used for photolysis.  
! This flexibility can be useful but users should be aware.
! For xsections that are temperature dependent, caution should be used in passing the proper 
! temperature to the data routines.  Usually, transmission is for layers, TLAY(NZ-1), while
! photolysis is at levels, T(NZ).
! The following subroutines are her: 
!     rdo3xs
!       o3_mol
!       o3_rei
!       o3_bas
!       o3_wmo
!       o3_jpl
!     rdo2xs
!     rdno2xs
!       no2xs_d
!       no2xs_jpl94
!       no2xs_har
!       no2xs_jpl06a
!       no2xs_jpl06b
!     rdso2xs
!=============================================================================*

      module module_xsections

      use phot_kind_mod, only: rk => kind_phot
      use params_mod, only: deltax, kin, input_data_root, qnan
      use  numer_mod, only: addpnt, inter2
      
      IMPLICIT NONE

      public :: o3xs,  no2xs_jpl06a
      public :: rdxs_init
      public :: o2_xs, so2_xs
      
      private

      REAL(rk), allocatable :: rei218(:), rei228(:), rei243(:), rei295(:)
      REAL(rk) :: v195, v345, v830
      REAL(rk), allocatable :: wmo203(:), wmo273(:)
      REAL(rk) :: v176, v850

      REAL(rk), allocatable :: jpl295(:), jpl218(:)
      REAL(rk) :: v186, v825

      REAL(rk), allocatable :: mol226(:), mol263(:), mol298(:)
      REAL(rk) :: v185, v240, v350

      REAL(rk), allocatable :: c0(:), c1(:), c2(:)
      REAL(rk) vb245, vb342

      REAL(rk), allocatable :: no2xs_a(:), no2xs_b(:)

      real(rk), protected, allocatable :: o2_xs(:)
      real(rk), protected, allocatable :: so2_xs(:)

    CONTAINS

      SUBROUTINE rdxs_init( nw, wl, errmsg, errflg )

      integer, intent(in) :: nw
      real(rk), intent(in) :: wl(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: istat, astat
      real(rk) :: xnan
      xnan = qnan()

      errmsg = ' '
      errflg = 0

      istat = 0
      if( .not. allocated( rei218 ) ) then
         allocate( rei218(nw),rei228(nw),rei243(nw),rei295(nw),stat=astat )
         istat = istat + astat
         rei218 = xnan; rei228=xnan; rei243=xnan; rei295=xnan
      endif
      if( .not. allocated( wmo203 ) ) then
         allocate( wmo203(nw),wmo273(nw),stat=astat )
         wmo203 = xnan; wmo273=xnan
         istat = istat + astat
      endif
      if( .not. allocated( jpl218 ) ) then
         allocate( jpl218(nw),jpl295(nw),stat=astat )
         jpl218=xnan; jpl295=xnan
         istat = istat + astat
      endif
      if( .not. allocated( mol226 ) ) then
         allocate( mol226(nw),mol263(nw),mol298(nw),stat=astat )
         mol226=xnan; mol263=xnan; mol298=xnan
         istat = istat + astat
      endif
      if( .not. allocated( c0 ) ) then
         allocate( c0(nw),c1(nw),c2(nw),stat=astat )
         c0=xnan; c1=xnan; c2=xnan
         istat = istat + astat
      endif
      if( .not. allocated( no2xs_a ) ) then
         allocate( no2xs_a(nw),no2xs_b(nw),stat=astat )
         no2xs_a=xnan; no2xs_b=xnan
         istat = istat + astat
      endif
      if (.not. allocated(o2_xs) ) then
         allocate(o2_xs(nw),stat=astat)
         o2_xs = xnan
         istat = istat + astat
      endif
      if (.not. allocated(so2_xs) ) then
         allocate(so2_xs(nw),stat=astat)
         o2_xs = xnan
         istat = istat + astat
      endif
      if( istat /= 0 ) then
         write(errmsg,'(''rdxs_init: failed to allocate; error = '',i4)') astat
         errflg = 1
         return
      endif
     
!_______________________________________________________________________
! read data from different sources
! rei = Reims group (Malicet et al., Brion et al.)
! jpl = JPL 2006 evaluation
! wmo = WMO 1985 O3 assessment
! mol = Molina and Molina
! bas = Bass et al.
!_______________________________________________________________________
      CALL o3_rei(nw,wl, errmsg, errflg)
      if (errflg.ne.0) return
      CALL o3_jpl(nw,wl, errmsg, errflg)
      if (errflg.ne.0) return
      CALL o3_wmo(nw,wl, errmsg, errflg)
      if (errflg.ne.0) return
      CALL o3_mol(nw,wl, errmsg, errflg)
      if (errflg.ne.0) return
      CALL o3_bas(nw,wl, errmsg, errflg)
      if (errflg.ne.0) return

      CALL rdno2xs(nw+1,wl, errmsg, errflg)
      call rdo2xs(nw+1,wl,o2_xs, errmsg, errflg)
      call rdso2xs(nw+1,wl,so2_xs, errmsg, errflg)
      
      END SUBROUTINE rdxs_init

      SUBROUTINE o3xs(nz,t,nw,wl, xs)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read ozone molecular absorption cross section.  Re-grid data to match    =*
!=  specified wavelength working grid. Interpolate in temperature as needed  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  MABS   - INTEGER, option for splicing different combinations of       (I)=*
!=           absorption cross secttions                                      =*
!=  NZ     - INTEGER, number of altitude levels or layers                 (I)=*
!=  T      - REAL, temperature of levels or layers                        (I)=*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid. In vacuum, nm                          =*
!=  XS     - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
!=           each specified wavelength (WMO value at 273)                    =*
!-----------------------------------------------------------------------------*

! input: (altitude working grid)

      INTEGER, intent(in) :: nw
      REAL(rk), intent(in)    :: wl(nw)

      INTEGER, intent(in) :: nz
      REAL(rk), intent(in) :: t(nz)

! output:
! ozone absorption cross sections interpolated to 
!   working wavelength grid (iw)
!   working altitude grid (iz) for temperature of layer or level (specified in call)
! Units are cm2 molecule-1 in vacuum

      REAL(rk), intent(inout) :: xs(:,:)

! internal

      INTEGER :: iw
      REAL(rk)    :: factor
      xs = 0._rk
      
!***** option 1:
! assign according to wavelength range:
!  175.439 - 185.185  1985WMO (203, 273 K)
!  185.185 - 195.00   2006JPL_O3 (218, 295 K)
!  195.00  - 345.00   Reims group (218, 228, 243, 295 K)
!  345.00  - 830.00   Reims group (295 K)
!  no extrapolations in temperature allowed

      DO iw = 1, nw-1
        IF(wl(iw) < v185) THEN
          factor = (wmo273(iw) - wmo203(iw))/(273._rk - 203._rk)
          xs(1:nz,iw) = wmo203(iw) + (t(1:nz) - 203._rk)*factor
          WHERE (t(1:nz) <= 203._rk) 
            xs(1:nz,iw) = wmo203(iw)
          ELSEWHERE (t(1:nz) >= 273._rk) 
            xs(1:nz,iw) = wmo273(iw)
          ENDWHERE
        ELSEIF(wl(iw) >= v185 .AND. wl(iw) < v195) THEN
          factor = (jpl295(iw) - jpl218(iw))/(295._rk - 218._rk)
          xs(1:nz,iw) = jpl218(iw) + (t(1:nz) - 218._rk)*factor
          WHERE (t(1:nz) <= 218._rk) 
            xs(1:nz,iw) = jpl218(iw)
          ELSEWHERE (t(1:nz) >= 295._rk) 
            xs(1:nz,iw) = jpl295(iw)
          ENDWHERE
        ELSEIF(wl(iw) >= v195 .AND. wl(iw) < v345) THEN
          factor = .1_rk*(rei228(iw) - rei218(iw))
          WHERE( t(1:nz) < 218._rk ) 
            xs(1:nz,iw) = rei218(iw)
          ELSEWHERE( t(1:nz) >= 218._rk .AND. t(1:nz) < 228._rk )
            xs(1:nz,iw) = rei218(iw) + (t(1:nz) - 218._rk)*factor
          ENDWHERE
          factor = (rei243(iw) - rei228(iw))/15._rk
          WHERE( t(1:nz) >= 228._rk .AND. t(1:nz) < 243._rk )
            xs(1:nz,iw) = rei228(iw) + (t(1:nz) - 228._rk)*factor
          ENDWHERE
          factor = (rei295(iw) - rei243(iw))/(295._rk - 243._rk)
          WHERE( t(1:nz) >= 243._rk .AND. t(1:nz) < 295._rk)
            xs(1:nz,iw) = rei243(iw) + (t(1:nz) - 243._rk)*factor
          ELSEWHERE( t(1:nz) >= 295._rk )
            xs(1:nz,iw) = rei295(iw)
          ENDWHERE
        ELSEIF(wl(iw) >= v345) THEN
          xs(1:nz,iw) = rei295(iw)
        ENDIF
      END DO

      END SUBROUTINE o3xs

!=============================================================================*

      SUBROUTINE o3_rei(nw,wl, errmsg, errflg)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and interpolate the O3 cross section from Reims group               =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  REI218 - REAL, cross section (cm^2) for O3 at 218K                    (O)=*
!=  REI228 - REAL, cross section (cm^2) for O3 at 218K                    (O)=*
!=  REI243 - REAL, cross section (cm^2) for O3 at 218K                    (O)=*
!=  REI295 - REAL, cross section (cm^2) for O3 at 218K                    (O)=*
!=  V195   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!=              e.g. start, stop, or other change                            =*
!=  V345   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!=  V830   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!-----------------------------------------------------------------------------*

!  input

      INTEGER, intent(in) :: nw
      REAL(rk), intent(in)    :: wl(nw)

! output:

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!* internal

      INTEGER, PARAMETER :: kdata = 70000

      INTEGER n1, n2, n3, n4
      REAL(rk) x1(kdata), x2(kdata), x3(kdata), x4(kdata)
      REAL(rk) y1(kdata), y2(kdata), y3(kdata), y4(kdata)

      INTEGER i
      INTEGER ierr

! used for air-to-vacuum wavelength conversion

      REAL(rk) ri(kdata)

      errmsg = ' '
      errflg = 0

! data from the Reims group:
!=  For Hartley and Huggins bands, use temperature-dependent values from     =*
!=  Malicet et al., J. Atmos. Chem.  v.21, pp.263-273, 1995.                 =*
!=  over 345.01 - 830.00, use values from Brion, room temperature only

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/O3/1995Malicet_O3.txt',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
         errmsg =  'o3_rei: Failed to open DATAE1/O3/1985Malicet_O3.txt'
         errflg = ierr
         return
      endif
      DO i = 1, 2
         READ(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           errmsg = 'o3_rei: Failed to read DATAE1/O3/1985Malicet_O3.txt'
           errflg = ierr
           return
         endif
      ENDDO
      n1 = 15001
      n2 = 15001
      n3 = 15001
      n4 = 15001
      DO i = 1, n1
         READ(kin,*,iostat=ierr) x1(i), y1(i), y2(i), y3(i), y4(i)
         if( ierr /= 0 ) then
           errmsg = 'o3_rei: Failed to read DATAE1/O3/1985Malicet_O3.txt' 
           errflg = ierr
           return
         endif
         x2(i) = x1(i)
         x3(i) = x1(i)
         x4(i) = x1(i)
      ENDDO
      CLOSE (kin)

!=  over 345.01 - 830.00, use values from Brion, room temperature only
! skip datum at 345.00 because already read in from 1995Malicet

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/O3/1998Brion_295.txt',STATUS='old')
      DO i = 1, 15
         READ(kin,*)
      ENDDO
      DO i = 1, 48515-15
         n1 = n1 + 1
         READ(kin,*) x1(n1), y1(n1)
      ENDDO
      CLOSE (kin)

      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19_rk)
      ENDDO
      DO i = 1, n1
         x1(i) = x1(i) * ri(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,rei295,n1,x1,y1,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_rei: interp err = '',i5,'' in O3 xsect - Reims 295K'')') ierr
         return
      ENDIF

      DO i = 1, n2
         ri(i) = refrac(x2(i), 2.45E19_rk)
      ENDDO
      DO i = 1, n2
         x2(i) = x2(i) * ri(i)
         x3(i) = x2(i)
         x4(i) = x2(i)
      ENDDO

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,rei243,n2,x2,y2,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_wmo: interp err = '',i5,'' in O3 xsect - Reims 243K'')') ierr
         return
      ENDIF

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,rei228,n3,x3,y3,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_wmo: interp err = '',i5,'' in O3 xswct - Reims 228K'')') ierr
         return
      ENDIF

      CALL addpnt(x4,y4,kdata,n4,x4(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x4,y4,kdata,n4,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x4,y4,kdata,n4,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,rei218,n4,x4,y4,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_wmo: interp err = '',i5,'' in O3 xswct - Reims 218K'')') ierr
         return
     ENDIF

! wavelength breaks must be converted to vacuum:

      v195 = 195.00_rk * refrac(195.00_rk, 2.45E19_rk)
      v345 = 345.00_rk * refrac(345.00_rk, 2.45E19_rk)
      v830 = 830.00_rk * refrac(830.00_rk, 2.45E19_rk)

      END SUBROUTINE o3_rei

!=============================================================================*

      SUBROUTINE o3_wmo(nw,wl, errmsg, errflg)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and interpolate the O3 cross section                                =*
!=  data from WMO 85 Ozone Assessment                                        =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WMO203 - REAL, cross section (cm^2) for O3 at 203K                    (O)=*
!=  WMO273 - REAL, cross section (cm^2) for O3 at 273K                    (O)=*
!=  V176   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!=              e.g. start, stop, or other change                            =*
!=  V850   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!-----------------------------------------------------------------------------*

!  input

      INTEGER, intent(in) :: nw
      REAL(rk), intent(in)    :: wl(nw)

! output
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! internal

      INTEGER, parameter :: kdata = 200

      INTEGER n1, n2
      REAL(rk) x1(kdata), x2(kdata)
      REAL(rk) y1(kdata), y2(kdata)

      INTEGER i, idum
      REAL(rk) a1, a2, dum
      INTEGER ierr

! used for air-to-vacuum wavelength conversion

      REAL(rk) ri(kdata)

! output

      errmsg = ' '
      errflg = 0
!----------------------------------------------------------
! cross sections from WMO 1985 Ozone Assessment
! from 175.439 to 847.500 nm

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/wmo85',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
         errmsg = 'o3_wmo: Failed to open DATAE1/wmo85'
         errflg = ierr
         return
      endif
      DO i = 1, 3
         read(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           errmsg = 'o3_wmo: Failed to open DATAE1/wmo85'
           errflg = ierr
           return
         endif
      ENDDO
      n1 = 158
      n2 = 158
      DO i = 1, n1
         READ(kin,*,iostat=ierr) idum, a1, a2, dum, dum, dum, y1(i), y2(i)
         if( ierr /= 0 ) then
           errmsg = 'o3_wmo: Failed to open DATAE1/wmo85'
           errflg = ierr
           return
         endif
         x1(i) = (a1+a2)/2._rk
         x2(i) = (a1+a2)/2._rk
      ENDDO
      CLOSE (kin)

! convert wavelengths to vacuum

      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19_rk)
      ENDDO
      DO i = 1, n1
         x1(i) = x1(i) * ri(i)
         x2(i) = x2(i) * ri(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,           1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,wmo203,n1,x1,y1,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_wmo: interp err = '',i5,'' in O3 cross section - WMO - 203K'')') ierr
         return
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,           1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,wmo273,n2,x2,y2,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_wmo: interp err = '',i5,'' in O3 cross section - WMO - 273K'')') ierr
         return
      ENDIF

! wavelength breaks must be converted to vacuum:
      
      a1 = (175.438_rk + 176.991_rk) / 2._rk
      v176 = a1 * refrac(a1,2.45E19_rk)

      a1 = (847.5_rk + 852.5_rk) / 2._rk
      v850 = a1 * refrac(a1, 2.45E19_rk)

      END SUBROUTINE o3_wmo

!=============================================================================*

      SUBROUTINE o3_jpl(nw,wl, errmsg, errflg)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and interpolate the O3 cross section from JPL 2006                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  JPL218 - REAL, cross section (cm^2) for O3 at 218K                    (O)=*
!=  JPL295 - REAL, cross section (cm^2) for O3 at 295K                    (O)=*
!=  V186   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!=              e.g. start, stop, or other change                            =*
!=  V825   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!-----------------------------------------------------------------------------*

!  input

      INTEGER, intent(in) :: nw
      REAL(rk), intent(in)    :: wl(nw)

! output:

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! internal

      INTEGER, parameter :: kdata = 200

      INTEGER n1, n2
      REAL(rk) x1(kdata), x2(kdata)
      REAL(rk) y1(kdata), y2(kdata)

      INTEGER i
      REAL(rk) dum
      INTEGER ierr

! used for air-to-vacuum wavelength conversion

      REAL(rk) ri(kdata)

      errmsg = ' '
      errflg = 0
      
! output

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/O3/2006JPL_O3.txt',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
         errmsg = 'o3_jpl: Failed to open DATAE1/O3/2006JPL_O3.txt'
         errflg = ierr
         return
      endif
      DO i = 1, 2
         read(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           errmsg = 'o3_jpl: Failed to read DATAE1/O3/2006JPL_O3.txt'
           errflg = ierr
           return
         endif
      ENDDO
      n1 = 167
      n2 = 167
      DO i = 1, n1
         READ(kin,*,iostat=ierr) dum, dum, x1(i), y1(i), y2(i)
         if( ierr /= 0 ) then
           errmsg = 'o3_jpl: Failed to read DATAE1/O3/2006JPL_O3.txt'
           errflg = ierr
           return
         endif
         y1(i) = y1(i) * 1.e-20_rk
         y2(i) = y2(i) * 1.e-20_rk
      ENDDO
      CLOSE (kin)

! convert wavelengths to vacuum

      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19_rk)
      ENDDO
      DO i = 1, n1
         x1(i) = x1(i) * ri(i)
         x2(i) = x1(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,           1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,jpl295,n1,x1,y1,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_jpl: interp err = '',i5,'' in file O3 cross section - WMO - 295K'')') ierr
         return
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,           1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,jpl218,n2,x2,y2,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_jpl: interp err = '',i5,'' in file O3 cross section - WMO - 218K'')') ierr
         return
      ENDIF

! wavelength breaks must be converted to vacuum:

      v186 = 186.051_rk * refrac(186.051_rk, 2.45E19_rk)
      v825 = 825._rk    * refrac(825._rk   , 2.45E19_rk)

      END SUBROUTINE o3_jpl

!=============================================================================*

      SUBROUTINE o3_mol(nw,wl, errmsg, errflg)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and interpolate the O3 cross section from Molina and Molina 1986    =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  MOL226 - REAL, cross section (cm^2) for O3 at 226 K                   (O)=*
!=  MOL263 - REAL, cross section (cm^2) for O3 at 263 K                   (O)=*
!=  MOL298 - REAL, cross section (cm^2) for O3 at 298 K                   (O)=*
!=  V185   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!=              e.g. start, stop, or other change                            =*
!=  V240   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!=  V350   - REAL, exact wavelength in vacuum for data breaks             (O)=*
!-----------------------------------------------------------------------------*

!  input

      INTEGER, intent(in) :: nw
      REAL(rk), intent(in)    :: wl(nw)

! output:

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! internal

      INTEGER i
      INTEGER ierr

      INTEGER, parameter :: kdata = 335
      INTEGER n1, n2, n3
      REAL(rk) x1(kdata), x2(kdata), x3(kdata)
      REAL(rk) y1(kdata), y2(kdata), y3(kdata)

! used for air-to-vacuum wavelength conversion

      REAL(rk) ri(kdata)

      errmsg = ' '
      errflg = 0

! output

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/O3/1986Molina.txt',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
         errmsg =  'o3_mol: Failed to open DATAE1/O3/1986Molina.txt'
         errflg = ierr
         return
      endif
      DO i = 1, 10
        READ(kin,*,iostat=ierr)
        if( ierr /= 0 ) then
          errmsg =  'o3_mol: Failed to read DATAE1/O3/1986Molina.txt'
          errflg = ierr
          return
       endif
      ENDDO
      n1 = 0
      n2 = 0
      n3 = 0
      DO i = 1, 121-10
         n1 = n1 + 1
         n3 = n3 + 1
         READ(kin,*,iostat=ierr) x1(n1), y1(n1),  y3(n3)
         if( ierr /= 0 ) then
            errmsg =  'o3_mol: Failed to read DATAE1/O3/1986Molina.txt'
            errflg = ierr
            return
         endif
         x3(n3) = x1(n1)
      ENDDO
      DO i = 1, 341-122
         n1 = n1 + 1
         n2 = n2 + 1
         n3 = n3 + 1
         READ(kin,*,iostat=ierr) x1(n1), y1(n1), y2(n2), y3(n3)
         if( ierr /= 0 ) then
            errmsg =  'o3_mol: Failed to read DATAE1/O3/1986Molina.txt'
            errflg = ierr
            return
         endif
         x2(n2) = x1(n1)
         x3(n3) = x1(n1)
      ENDDO
      CLOSE (kin)

! convert all wavelengths from air to vacuum

      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19_rk)
      ENDDO
      DO i = 1, n1
         x1(i) = x1(i) * ri(i)
      ENDDO

      DO i = 1, n2
         ri(i) = refrac(x2(i), 2.45E19_rk)
      ENDDO
      DO i = 1, n2
         x2(i) = x2(i) * ri(i)
      ENDDO

      DO i = 1, n3
         ri(i) = refrac(x3(i), 2.45E19_rk)
      ENDDO
      DO i = 1, n3
         x3(i) = x3(i) * ri(i)
      ENDDO

! convert wavelength breaks from air to vacuum

      v185 = 185._rk  * refrac(185._rk , 2.45E19_rk)
      v240 = 240.5_rk * refrac(240.5_rk, 2.45E19_rk)
      v350 = 350._rk  * refrac(350._rk , 2.45E19_rk)

! interpolate to working grid

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,mol226,n1,x1,y1,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_mol: interp err = '',i5,'' in O3 xsect - 226K Molina'')') ierr
         return
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,mol263,n2,x2,y2,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_mol: interp err = '',i5,'' in O3 xsect - 263K Molina'')') ierr
         return
      ENDIF

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,mol298,n3,x3,y3,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_mol: interp err = '',i5,'' in O3 xsect - 298K Molina'')') ierr
         return
      ENDIF

      END SUBROUTINE o3_mol

!=============================================================================*

      SUBROUTINE o3_bas(nw,wl, errmsg, errflg)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and interpolate the O3 cross section from Bass 1985                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  c0     - REAL, coefficint for polynomial fit to cross section (cm^2)  (O)=*
!=  c1     - REAL, coefficint for polynomial fit to cross section (cm^2)  (O)=*
!=  c2     - REAL, coefficint for polynomial fit to cross section (cm^2)  (O)=*
!=  Vb245   - REAL, exact wavelength in vacuum for data breaks            (O)=*
!=              e.g. start, stop, or other change                            =*
!=  Vb342   - REAL, exact wavelength in vacuum for data breaks            (O)=*
!-----------------------------------------------------------------------------*

! input:

      INTEGER, intent(in) :: nw
      REAL(rk), intent(in)    :: wl(nw)

! output:

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! internal:

      INTEGER, parameter :: kdata = 2000

      INTEGER i
      INTEGER ierr

      INTEGER n1, n2, n3
      REAL(rk) x1(kdata), x2(kdata), x3(kdata)
      REAL(rk) y1(kdata), y2(kdata), y3(kdata)

! used for air-to-vacuum wavelength conversion

      REAL(rk) ri(kdata)

      errmsg = ' '
      errflg = 0

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/O3/1985Bass_O3.txt',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
          errmsg = 'o3_bas: Failed to open DATAE1/O3/1985Bass_O3.txt'
          errflg = ierr
          return
      endif
      DO i = 1, 8
         READ(kin,*,iostat=ierr)
      ENDDO
      n1 = 1915
      n2 = 1915
      n3 = 1915
      DO i = 1, n1
        READ(kin,*,iostat=ierr) x1(i), y1(i), y2(i), y3(i)
        if( ierr /= 0 ) then
          errmsg = 'o3_bas: Failed to read DATAE1/O3/1985Bass_O3.txt'
          errflg = ierr
          return
        endif
      ENDDO
      CLOSE (kin)
      y1(1:n1) = 1.e-20_rk * y1(1:n1)
      y2(1:n1) = 1.e-20_rk * y2(1:n1)
      y3(1:n1) = 1.e-20_rk * y3(1:n1)

! convert all wavelengths from air to vacuum

      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19_rk)
      ENDDO
      x1(1:n1) = x1(1:n1) * ri(1:n1)
      x2(1:n1) = x1(1:n1)
      x3(1:n1) = x1(1:n1)

! convert wavelength breaks to vacuum

      vb245 = 245.018_rk * refrac(245.018_rk, 2.45E19_rk)
      vb342 = 341.981_rk * refrac(341.981_rk, 2.45E19_rk)

! interpolate to working grid

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,c0,n1,x1,y1,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_bas: interp err = '',i5,'' in O3 xsect - c0 Bass'')') ierr
         return
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,c1,n2,x2,y2,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_bas: interp err = '',i5,'' in O3 xsect - c1 Bass'')') ierr
         return
      ENDIF

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x3,y3,kdata,n3,            1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,c2,n3,x3,y3,errmsg, errflg)
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''o3_bas: interp err = '',i5,'' in O3 xsect - c2 Bass'')') ierr
         return
      ENDIF

      END SUBROUTINE o3_bas

!=============================================================================*

      SUBROUTINE rdo2xs(nw,wl,o2xs1, errmsg, errflg)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute equivalent O2 cross section, except                              =*
!=  the SR bands and the Lyman-alpha line.                                   =*
!-----------------------------------------------------------------------------* 
!=  PARAMETERS:                                   
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
!=            wavelength grid                                                =*
!=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
!=            working wavelength grid           
!=            vertical layer at each specified wavelength                    =*
!=  O2XS1   - REAL, O2 molecular absorption cross section                    =*
!=
!-----------------------------------------------------------------------------*

! Input

      INTEGER, intent(in)  :: nw
      REAL(rk), intent(in) :: wl(:)

! Output O2 xsect, temporary, will be over-written in Lyman-alpha and 
!   Schumann-Runge wavelength bands.

      REAL(rk), intent(inout) :: o2xs1(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! Internal

      INTEGER, parameter :: kdata = 200
      INTEGER :: i, n
      INTEGER :: ierr
      REAL(rk)    :: x, y
      REAL(rk)    :: x1(kdata), y1(kdata)

      errmsg = ' '
      errflg = 0

! Read O2 absorption cross section data:
!  116.65 to 203.05 nm = from Brasseur and Solomon 1986
!  205 to 240 nm = Yoshino et al. 1988

! Note that subroutine la_srb.f will over-write values in the spectral regions
!   corresponding to:
! - Lyman-alpha (LA: 121.4-121.9 nm, Chabrillat and Kockarts parameterization) 
! - Schumann-Runge bands (SRB: 174.4-205.8 nm, Koppers parameteriaztion)

      n = 0

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/O2/O2_brasseur.abs',iostat=ierr)
      if( ierr /= 0 ) then
         errmsg = 'rdo2xs: Failed to open DATAE1/O2/O2_brasseur.abs'
         errflg = ierr
         return
      endif
      DO i = 1, 7
         READ(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           errmsg = 'rdo2xs: Failed to open DATAE1/O2/O2_brasseur.abs'
           errflg = ierr
           return
         endif
      ENDDO
      DO i = 1, 78
         READ(kin,*,iostat=ierr) x, y
         if( ierr /= 0 ) then
           errmsg = 'rdo2xs: Failed to open DATAE1/O2/O2_brasseur.abs'
           errflg = ierr
           return
         endif
         IF (x .LE. 204._rk) THEN
            n = n + 1
            x1(n) = x
            y1(n) = y
         ENDIF
      ENDDO
      CLOSE(kin)

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/O2/O2_yoshino.abs',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
         errmsg = 'rdo2xs: Failed to open DATAE1/O2/O2_brasseur.abs'
         errflg = ierr
         return
      endif

      DO i = 1, 8
         READ(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
            errmsg = 'rdo2xs: Failed to read DATAE1/O2/O2_yoshino.abs'
            errflg = ierr
            return
         endif
      ENDDO
      DO i = 1, 36
         n = n + 1
         READ(kin,*,iostat=ierr) x, y
         if( ierr /= 0 ) then
            errmsg = 'rdo2xs: Failed to read DATAE1/O2/O2_yoshino.abs'
            errflg = ierr
            return
         endif
         y1(n) = y*1.E-24_rk
         x1(n) = x
      END DO
      CLOSE (kin)

! Add termination points and interpolate onto the 
!  user grid (set in subroutine gridw):

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1._rk-deltax),y1(1),errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n,0._rk               ,y1(1),errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n,              1.E+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,o2xs1, n,x1,y1,errmsg, errflg)
      
      IF (errflg .NE. 0) THEN
         WRITE(errmsg,'(''rdo2xs: interp err = '',i5,'' in O2 -> O + O'')') ierr
         return
      ENDIF

      END SUBROUTINE rdo2xs

!=============================================================================*

      SUBROUTINE rdno2xs(nw,wl, errmsg, errflg)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read NO2 molecular absorption cross section.  Re-grid data to match      =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  NO2XS  - REAL, molecular absoprtion cross section (cm^2) of NO2 at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*

! input:

      INTEGER, intent(in)  :: nw
      REAL(rk), intent(in) :: wl(:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! locals:
      INTEGER, parameter :: kdata = 100
      INTEGER :: i, n1, n2
      REAL(rk)    :: dum1, dum2
      REAL(rk)    :: x1(kdata), x2(kdata), y1(kdata), y2(kdata)

      errmsg = ' '
      errflg = 0

! NO2 absorption cross section from JPL2006
! with interpolation of bin midpoints

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/NO2/NO2_jpl2006.abs',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO 
      n1 = 81
      DO i = 1, n1
         READ(kin,*) dum1, dum2, y1(i), y2(i)
         x1(i) = 0.5_rk * (dum1 + dum2)
         x2(i) = x1(i) 
         y1(i) = y1(i)*1.E-20_rk
         y2(i) = y2(i)*1.E-20_rk
      ENDDO
      CLOSE(kin)
      n2 = n1

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1._rk+deltax),   0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n1,            1.e+38_rk,   0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,no2xs_a,n1,x1,y1,errmsg, errflg)
      if ( errflg .ne. 0) return
      
      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,               0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1._rk+deltax),   0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x2,y2,kdata,n2,            1.e+38_rk,   0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,no2xs_b,n2,x2,y2,errmsg, errflg)
      if ( errflg .ne. 0) return

      END SUBROUTINE rdno2xs

!=============================================================================*

      SUBROUTINE no2xs_jpl06a(nz,t,nw,wl, no2xs)

! interpolate NO2 xs from JPL2006

! input:

      INTEGER, intent(in) :: nz
      INTEGER, intent(in) :: nw
      REAL(rk), intent(in)    :: t(nz)
      REAL(rk), intent(in)    :: wl(nw)

! output:

      REAL(rk), intent(inout) :: no2xs(:,:)

! local

      INTEGER :: iw
      REAL(rk)    :: tfac(nz)
      
      tfac(1:nz) = (t(1:nz) - 220._rk)/74._rk
      DO iw = 1, nw-1
        no2xs(1:nz,iw) = no2xs_a(iw) + (no2xs_b(iw)-no2xs_a(iw))*tfac(1:nz)
      ENDDO 

      END SUBROUTINE no2xs_jpl06a

!=============================================================================*

      SUBROUTINE rdso2xs(nw,wl,so2xs, errmsg, errflg)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read SO2 molecular absorption cross section.  Re-grid data to match      =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*

! input: (altitude working grid)
      INTEGER, intent(in)  :: nw
      REAL(rk), intent(in) :: wl(:)

! output:

      REAL(rk), intent(inout) :: so2xs(:)

      character(len=*), intent(out)   :: errmsg
      integer,          intent(out)   :: errflg

!! local:
      INTEGER, parameter :: kdata = 1000

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)
      INTEGER i, n
      CHARACTER(len=40)  :: fil

      
      errmsg = ' '
      errflg = 0
!************ absorption cross sections:
! SO2 absorption cross sections from J. Quant. Spectrosc. Radiat. Transfer
! 37, 165-182, 1987, T. J. McGee and J. Burris Jr.
! Angstrom vs. cm2/molecule, value at 221 K

      fil = 'DATA/McGee87'
      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAE1/SO2/SO2xs.all',STATUS='old')
      DO i = 1,3 
        read(kin,*)
      ENDDO
      n = 704 
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        x1(i) = .1_rk*x1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n,          0._rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1._rk+deltax),0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL addpnt(x1,y1,kdata,n,      1.e+38_rk,0._rk,errmsg, errflg)
      if (errflg.ne.0) return
      CALL inter2(nw,wl,so2xs,n,x1,y1,errmsg, errflg)

      END SUBROUTINE rdso2xs

      real(rk) FUNCTION refrac(w,airden)

      IMPLICIT NONE

! input vacuum wavelength, nm and air density, molec cm-3

      REAL(rk), intent(in) :: w, airden

! output refractive index for standard air
! (dry air at 15 deg. C, 101.325 kPa, 0.03% CO2)

! internal

      REAL(rk) :: sig,  sigsq, dum

! from CRC Handbook, originally from Edlen, B., Metrologia, 2, 71, 1966.
! valid from 200 nm to 2000 nm
! beyond this range, use constant value

      IF (w < 200._rk) then
        dum = 5.e-3_rk
      ELSEIF (w > 2000._rk) then
        dum = 5.e-4_rk
      ELSE
        dum = 1._rk/w
      ENDIF
      sig = 1.E3_rk*dum
      sigsq = sig * sig

      dum = 8342.13_rk + 2406030._rk/(130._rk - sigsq) + 15997._rk/(38.9_rk - sigsq)

! adjust to local air density
      dum = dum * airden/(2.69e19_rk * 273.15_rk/288.15_rk)

! index of refraction:
      refrac = 1._rk + 1.E-8_rk * dum

      END function refrac

      end module module_xsections
