module phot_util_mod
  use phot_kind_mod, only: dp
  use phot_kind_mod, only: rk => kind_phot

  implicit none
  public
contains
  
      REAL(rk) FUNCTION sundis( julday )
!-----------------------------------------------------------------------------
! purpose:
! calculate earth-sun distance variation for a given date. based on
! fourier coefficients originally from: spencer, j.w., 1971, fourier
! series representation of the position of the sun, search, 2:172
!-----------------------------------------------------------------------------
! parameters:
! idate - integer, specification of the date, from yymmdd (i)
! esrm2 - real(dp), variation of the earth-sun distance (o)
! esrm2 = (average e/s dist)^2 / (e/s dist on day idate)^2
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! ... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: julday

!-----------------------------------------------------------------------------
! ... local variables
!-----------------------------------------------------------------------------
      real(dp), parameter :: pi = 3.1415926_dp

      real(dp) :: dayn, thet0
      real(dp) :: sinth, costh, sin2th, cos2th

!-----------------------------------------------------------------------------
! ... parse date to find day number (julian day)
!-----------------------------------------------------------------------------
      dayn = real(julday - 1,kind=dp) + .5_dp
!-----------------------------------------------------------------------------
! ... define angular day number and compute esrm2:
!-----------------------------------------------------------------------------
      thet0 = 2._dp*pi*dayn/365._dp
!-----------------------------------------------------------------------------
! ... calculate sin(2*thet0), cos(2*thet0) from
! addition theorems for trig functions for better
! performance; the computation of sin2th, cos2th
! is about 5-6 times faster than the evaluation
! of the intrinsic functions sin and cos
!-----------------------------------------------------------------------------
      sinth  = sin( thet0 )
      costh  = cos( thet0 )
      sin2th = 2._dp*sinth*costh
      cos2th = costh*costh - sinth*sinth
      sundis  = real( 1.000110_dp + .034221_dp*costh + .001280_dp*sinth &
                                  + .000719_dp*cos2th + .000077_dp*sin2th )

      END FUNCTION sundis

      subroutine calc_zenith( lat, long, julday, gmt, zenith, &
                              its, ite, jts, jte, &
                              ims, ime, jms, jme )
!-------------------------------------------------------------------
! this subroutine calculates solar zenith and azimuth angles for a
! part  time and location.  must specify:
! input:
! lat - latitude in decimal degrees
! long - longitude in decimal degrees
! gmt  - greenwich mean time - decimal military eg.
! 22.75 = 45 min after ten pm gmt
! output:
! zenith
! azimuth
! .. Scalar Arguments ..
!-------------------------------------------------------------------
        integer,  intent(in)  :: julday
        integer,  intent(in)  :: its,ite
        integer,  intent(in)  :: jts,jte
        integer,  intent(in)  :: ims,ime
        integer,  intent(in)  :: jms,jme
        real(dp), intent(in)  :: gmt
        real(rk), intent(in)  :: lat(ims:ime,jms:jme)
        real(rk), intent(in)  :: long(ims:ime,jms:jme)
        real(rk), intent(out) :: zenith(ims:ime,jms:jme)

!-------------------------------------------------------------------
! .. Local variables
!-------------------------------------------------------------------
      real(dp), parameter :: d2r = 3.1415926_dp/180.0_dp
      real(dp), parameter :: r2d = 1.0_dp/d2r

      integer  :: i, j
      real(dp) :: csz, cw, d, ec, epsi, eqt, eyt, feqt, feqt1, &
          feqt2, feqt3, feqt4, feqt5, feqt6, feqt7, lbgmt, lzgmt, ml, pepsi, &
          ra,  rdecl, reqt, rlt, rml, rra, ssw, sw, tab, w, wr, &
          yt, zpt, zr

      d = real(julday,dp) + gmt/24.0_dp
!-------------------------------------------------------------------
! calc geom mean longitude
!-------------------------------------------------------------------
      ml  = 279.2801988_dp + d*(.9856473354_dp + 2.267E-13_dp*d)
      rml = ml*d2r
!-------------------------------------------------------------------
! calc equation of time in sec
! w = mean long of perigee
! e = eccentricity
! epsi = mean obliquity of ecliptic
!-------------------------------------------------------------------
      w = 282.4932328_dp + d*(4.70684E-5_dp + 3.39E-13_dp*d)
      wr = w*d2r
      ec   = 1.6720041E-2_dp - d*(1.1444E-9_dp + 9.4E-17_dp*d)
      epsi = 23.44266511_dp - d*(3.5626E-7_dp + 1.23E-15_dp*d)
      pepsi = epsi*d2r
      yt = (tan(pepsi/2.0_dp))**2
      cw = cos(wr)
      sw = sin(wr)
      ssw = sin(2.0_dp*wr)
      eyt = 2._dp*ec*yt
      feqt1 = -sin(rml)*cw*(eyt + 2._dp*ec)
      feqt2 = cos(rml)*sw*(2._dp*ec - eyt)
      feqt3 = sin(2._dp*rml)*(yt - (5._dp*ec**2/4._dp)*(cw**2 - sw**2))
      feqt4 = cos(2._dp*rml)*(5._dp*ec**2*ssw/4._dp)
      feqt5 = sin(3._dp*rml)*(eyt*cw)
      feqt6 = -cos(3._dp*rml)*(eyt*sw)
      feqt7 = -sin(4._dp*rml)*(.5_dp*yt**2)
      feqt  = feqt1 + feqt2 + feqt3 + feqt4 + feqt5 + feqt6 + feqt7
      eqt   = feqt*13751.0_dp

!-------------------------------------------------------------------
! convert eq of time from sec to deg
!-------------------------------------------------------------------
      reqt = eqt/240._dp
!-------------------------------------------------------------------
! calc right ascension in rads
!-------------------------------------------------------------------
      ra  = ml - reqt
      rra = ra*d2r
!-------------------------------------------------------------------
! calc declination in rads, deg
!-------------------------------------------------------------------
      tab = 0.43360_dp*sin(rra)
      rdecl = atan(tab)
      do j = jts,jte
        do i = its,ite
!-------------------------------------------------------------------
! calc local hour angle
!-------------------------------------------------------------------
          lbgmt = 12.0_dp - eqt/3600._dp + real(long(i,j),dp)*24._dp/360._dp
          lzgmt = 15.0_dp*(gmt - lbgmt)
          zpt   = lzgmt*d2r
          rlt   = real(lat(i,j),dp)*d2r
          csz   = sin(rlt)*sin(rdecl) + cos(rlt)*cos(rdecl)*cos(zpt)
          csz   = min( 1._dp,csz )
          zr    = acos(csz)
          zenith(i,j) = real( zr/d2r )
        end do
      end do

      end subroutine calc_zenith

      SUBROUTINE sphers( nlyr, z, zen, dsdh, nid )
!-----------------------------------------------------------------------------
!=  PURPOSE:
!=  Calculate slant path over vertical depth ds/dh in spherical geometry.
!=  Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model
!=  for computing the radiation field available for photolysis and heating
!=  at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)
!=            grid
!=  Z       - REAL, specified altitude working grid (km)                  (I)
!=  ZEN     - REAL, solar zenith angle (degrees)                          (I)
!=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)
!=            when travelling from the top of the atmosphere to layer i;
!=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1
!=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)
!=            travelling from the top of the atmosphere to layer i;
!=            NID(i), i = 0..NZ-1
!-----------------------------------------------------------------------------
!  This program is free software;  you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by the
!= Free Software Foundation;  either version 2 of the license, or (at your
!= option) any later version.
!= The TUV package is distributed in the hope that it will be useful, but
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
!= License for more details.
!= To obtain a copy of the GNU General Public License, write to:
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!-----------------------------------------------------------------------------
!= To contact the authors, please mail to:
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA
!= send email to:  sasha@ucar.edu
!-----------------------------------------------------------------------------

      use params_mod, only : pi, radius

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      INTEGER,  intent(in)    :: nlyr
      REAL(rk), intent(in)    :: zen
      REAL(rk), intent(in)    :: z(:)
      INTEGER,  intent(inout) :: nid(0:nlyr)
      REAL(rk), intent(inout) :: dsdh(0:nlyr,nlyr)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      REAL(rk), parameter ::  dr = pi/180._rk

      INTEGER :: j, jm1, k
      INTEGER :: id
      REAL(rk)    :: re
      REAL(rk)    :: zd(0:nlyr)
      REAL(rk)    :: ds_dh(1:nlyr)
      REAL(dp) :: zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm

      zenrad = REAL( zen*dr, dp )
!-----------------------------------------------------------------------------
! include the elevation above sea level to the radius of the earth:
!-----------------------------------------------------------------------------
      re = radius + z(1)

!-----------------------------------------------------------------------------
! from bottom-up to top-down
! note zd is the elevation above earth surface:
!-----------------------------------------------------------------------------
      zd(0:nlyr) = z(nlyr+1:1:-1) - z(1)

!-----------------------------------------------------------------------------
! initialize nid
!-----------------------------------------------------------------------------
      nid(0:nlyr) = 0

!-----------------------------------------------------------------------------
! calculate ds/dh of every layer
!-----------------------------------------------------------------------------
layer_loop : &
      DO k = 0, nlyr
        ds_dh(:) = 0._rk
        rpsinz = real(re + zd(k),dp) * SIN(zenrad)
!       IF( zen > 90.0 .AND. rpsinz < real(re,8) ) THEN
        IF( zen <= 90.0_rk .or. rpsinz >= real(re,dp) ) THEN
!-----------------------------------------------------------------------------
! Find index of layer in which the screening height lies
!-----------------------------------------------------------------------------
          id = k 
          IF( zen > 90.0_rk ) THEN
            DO j = 1,nlyr
              IF( rpsinz < real(zd(j-1) + re,dp) .AND. &
                  rpsinz >= real(zd(j) + re,dp) ) then
                id = j
              ENDIF
            END DO
          END IF
 
          DO j = 1, id
            jm1 = j - 1
!           IF( j == id .AND. id == k .AND. zen > 90.0 ) then
            IF( j /= id .or. k /= id .or. zen <= 90.0_rk ) then
              sm = 1.0_dp
            ELSE
              sm = -1.0_dp
            ENDIF
            rj   = real(re + zd(jm1),dp)
            rjp1 = real(re + zd(j),dp)
            dhj  = zd(jm1) - zd(j)
 
            ga = max( rj*rj - rpsinz*rpsinz,0.0_dp )
            gb = max( rjp1*rjp1 - rpsinz*rpsinz,0.0_dp )

            IF( id > k .AND. j == id ) THEN
              dsj = SQRT( ga )
            ELSE
              dsj = SQRT( ga ) - sm*SQRT( gb )
            END IF
            ds_dh(j) = real( dsj/dhj )
          END DO
          nid(k) = id
        ELSE
          nid(k) = -1
        ENDIF
        dsdh(k,:) = ds_dh(:)

      END DO layer_loop

      END SUBROUTINE sphers

      SUBROUTINE airmas( nlyr, dsdh, nid, cz, vcol, scol )
!-----------------------------------------------------------------------------
!=  PURPOSE:
!=  Calculate vertical and slant air columns, in spherical geometry, as a
!=  function of altitude.
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)
!=            grid
!=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)
!=            when travelling from the top of the atmosphere to layer i;
!=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1
!=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)
!=            travelling from the top of the atmosphere to layer i;
!=            NID(i), i = 0..NZ-1
!=  VCOL    - REAL, output, vertical air column, molec cm-2, above level iz
!=  SCOL    - REAL, output, slant air column in direction of sun, above iz
!=            also in molec cm-2
!-----------------------------------------------------------------------------

      use params_mod, only : largest

      IMPLICIT NONE

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      INTEGER, intent(in) :: nlyr
      INTEGER, intent(in) :: nid(0:nlyr)
      REAL(rk),    intent(in) :: dsdh(0:nlyr,nlyr)
      REAL(rk),    intent(in) :: cz(nlyr)

      REAL(rk), intent(inout) :: vcol(nlyr), scol(nlyr)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      INTEGER :: lyr, j, nlev, nlevi
      REAL(rk)    :: sum, vsum

!-----------------------------------------------------------------------------
! calculate vertical and slant column from each level: work downward
!-----------------------------------------------------------------------------
      nlev = nlyr + 1
      vsum = 0._rk
      DO lyr = 1, nlyr
        nlevi = nlev - lyr
        vsum = vsum + cz(nlevi)
        vcol(nlevi) = vsum
        sum = 0._rk
        IF( nid(lyr) < 0 ) THEN
          sum = largest
        ELSE
!-----------------------------------------------------------------------------
! single pass layers:
!-----------------------------------------------------------------------------
          DO j = 1, MIN(nid(lyr), lyr)
            sum = sum + cz(nlev-j)*dsdh(lyr,j)
          END DO
!-----------------------------------------------------------------------------
! double pass layers:
!-----------------------------------------------------------------------------
           DO j = MIN(nid(lyr),lyr)+1, nid(lyr)
             sum = sum + 2._rk*cz(nlev-j)*dsdh(lyr,j)
           END DO
        ENDIF
        scol(nlevi) = sum
      END DO

      END SUBROUTINE airmas

end module phot_util_mod
