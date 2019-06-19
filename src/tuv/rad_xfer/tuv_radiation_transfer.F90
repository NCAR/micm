module tuv_radiation_transfer
  use phot_kind_mod, only: rk => kind_phot

  implicit none

  private
  public :: tuv_radiation_transfer_init
  public :: tuv_radiation_transfer_run
  public :: tuv_radiation_transfer_finalize
  
  integer :: nlev, nlyr
  
contains
  
!> \section arg_table_tuv_radiation_transfer_init Argument Table
!! | local_name | standard_name             | long_name                 | units   | rank | type      | kind      | intent | optional |
!! |------------|---------------------------|---------------------------|---------|------|-----------|-----------|--------|----------|
!! | realkind   | phys_real_kind            | physics real kind         | none    |    0 | integer   |           | in     | F        |
!! | nlevels    | num_levels_for_photolysis | number of column layers   | count   |    0 | integer   |           | in     | F        |
!! | errmsg     | ccpp_error_message        | CCPP error message        | none    |    0 | character | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag           | CCPP error flag           | flag    |    0 | integer   |           | out    | F        |
!!
subroutine tuv_radiation_transfer_init( realkind, nlevels, errmsg, errflg )
    use params_mod,       only: input_data_root
    use rad_abs_xsect,    only: rad_abs_xsect_init
    use module_xsections, only: rdxs_init
    use rad_abs_xsect,    only: nwave, wl

    integer,          intent(in)  :: realkind
    integer,          intent(in)  :: nlevels
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=256) :: filepath
    
    errmsg = ''
    errflg = 0

    if ( realkind/=rk ) then
       errmsg = 'tuv_radiation_transfer_init: realkind does not match kind_phot'
       errflg = 1
       return
    end if

    nlev = nlevels
    nlyr = nlev-1

    filepath = trim(input_data_root)//'/wrf_tuv_xsqy.nc'

    call rad_abs_xsect_init( filepath, errmsg, errflg )
    if (errflg.ne.0) return
    
    call rdxs_init( nwave, wl, errmsg, errflg )
    if (errflg.ne.0) return
   
  end subroutine tuv_radiation_transfer_init
  
!> \section arg_table_tuv_radiation_transfer_run Argument Table
!! | local_name | standard_name                         | long_name                          | units     | rank | type      | kind      | intent | optional |
!! |------------|---------------------------------------|------------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | zenith     | solar_zenith                          | solar zenith angle                 | degrees   |    0 | real      | kind_phys | in     | F        |
!! | albedo     | surface_albedo                        | surface albedo                     | none      |    0 | real      | kind_phys | in     | F        |
!! | press_mid  | layer_pressure                        | mid-point layer pressure           | Pa        |    1 | real      | kind_phys | in     | F        |
!! | alt        | layer_altitude                        | mid-point layer altitude           | m         |    1 | real      | kind_phys | in     | F        |
!! | temp       | layer_temperature                     | mid-point layer temperature        | K         |    1 | real      | kind_phys | in     | F        |
!! | o3vmr      | O3_vmr_col                            | O3 volume mixing ratio column      | mole/mole |    1 | real      | kind_phys | in     | F        |
!! | so2vmr     | SO2_vmr_col                           | SO2 volume mixing ratio column     | mole/mole |    1 | real      | kind_phys | in     | F        |
!! | no2vmr     | NO2_vmr_col                           | NO2 volume mixing ratio column     | mole/mole |    1 | real      | kind_phys | in     | F        |
!! | dto2       | O2_optical_depth                      | optical depth due to O2 absorption | cm        |    2 | real      | kind_phys | in     | F        |
!! | radfld     | actinic_photon_fluxes                 | actinic photon fluxes              | cm-2 sec-1|    2 | real      | kind_phys | out    | F        |
!! | errmsg     | ccpp_error_message                    | CCPP error message                 | none      |    0 | character | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag                       | CCPP error flag                    | flag      |    0 | integer   |           | out    | F        |
!!
  subroutine tuv_radiation_transfer_run( zenith, albedo, press_mid, alt, temp, o3vmr, so2vmr, no2vmr, dto2, radfld, errmsg, errflg )

    use tuv_subs,         only: tuv_radfld
    use rad_abs_xsect,    only: o2_xs, so2_xs, nwave, wl, wc
    use module_xsections, only: o3xs, no2xs_jpl06a
    use params_mod
 
    real(rk),         intent(in)  :: zenith
    real(rk),         intent(in)  :: albedo
    real(rk),         intent(in)  :: press_mid(:)
    real(rk),         intent(in)  :: alt(:)  ! m
    real(rk),         intent(in)  :: temp(:) ! K
    real(rk),         intent(in)  :: o3vmr(:)
    real(rk),         intent(in)  :: so2vmr(:)
    real(rk),         intent(in)  :: no2vmr(:)
    real(rk),         intent(in)  :: dto2(:,:)
    real(rk),         intent(out) :: radfld(:,:) ! /sec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: nlambda_start=1
    integer, parameter :: cld_od_opt=1
    logical, parameter :: has_aer_ra_feedback = .false.
    real(rk), parameter :: dobsi = 0._rk
    
    real(rk) :: zen
    real(rk) :: alb(nwave)
    real(rk) :: zlev(nlev) ! km 
    real(rk) :: tlev(nlev)
    real(rk) :: aircol(nlyr)  ! # molecules / cm2 in each layer
    real(rk) :: o3col(nlyr) 
    real(rk) :: so2col(nlyr)
    real(rk) :: no2col(nlyr)
    real(rk) :: dpress(nlyr)

    real(rk) :: tauaer300(nlev) ! aerosol properties
    real(rk) :: tauaer400(nlev)
    real(rk) :: tauaer600(nlev)
    real(rk) :: tauaer999(nlev)
    real(rk) :: waer300(nlev)
    real(rk) :: waer400(nlev)
    real(rk) :: waer600(nlev)
    real(rk) :: waer999(nlev)
    real(rk) :: gaer300(nlev)
    real(rk) :: gaer400(nlev)
    real(rk) :: gaer600(nlev)
    real(rk) :: gaer999(nlev)
    
    real(rk) :: dtaer(nlyr,nwave), omaer(nlyr,nwave), gaer(nlyr,nwave)
    real(rk) :: dtcld(nlyr,nwave), omcld(nlyr,nwave), gcld(nlyr,nwave)
    real(rk) :: dt_cld(nlyr)
    
    real(rk) :: qll(nlev) ! cld water content (g/m3)
    real(rk) :: cldfrac(nlev)
    real(rk) :: efld(nlev,nwave)
    real(rk) :: e_dir(nlev,nwave)
    real(rk) :: e_dn(nlev,nwave)
    real(rk) :: e_up(nlev,nwave)
    real(rk) :: dir_fld(nlev,nwave)
    real(rk) :: dwn_fld(nlev,nwave)
    real(rk) :: up_fld(nlev,nwave)

    integer :: k, kk

    real(rk) :: o3_xs(nwave,nlev)
    real(rk) :: no2_xs(nwave,nlev)
    real(rk) :: o3_xs_tpose(nlev,nwave)
    real(rk) :: no2_xs_tpose(nlev,nwave)

    errmsg = ''
    errflg = 0

    dpress(1:nlyr) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    do k=1,nlyr
       kk=nlyr-k+1
       aircol(k) = 10._rk*dpress(k)*R/(kboltz*g)
       o3col(kk)  = 0.5_rk*(o3vmr(k)+o3vmr(k+1))*aircol(k)
       so2col(kk) = 0.5_rk*(so2vmr(k)+so2vmr(k+1))*aircol(k)
       no2col(kk) = 0.5_rk*(no2vmr(k)+no2vmr(k+1))*aircol(k)
    end do

    ! inputs need to be bottom up vert coord
    aircol(1:nlyr) = aircol(nlyr:1:-1)
    tlev(nlev:1:-1) = temp(1:nlev)
    zlev(nlev:1:-1) = alt(1:nlev)*1.e-3_rk ! m -> km

    qll=0.0_rk
    cldfrac=0.0_rk
    tauaer300=0.0_rk
    tauaer400=0.0_rk
    tauaer600=0.0_rk
    tauaer999=0.0_rk
    waer300=1.0_rk
    waer400=1.0_rk
    waer600=1.0_rk
    waer999=1.0_rk
    gaer300=0.0_rk
    gaer400=0.0_rk
    gaer600=0.0_rk
    gaer999=0.0_rk

    zen = zenith
    alb(:) = albedo

    o3_xs_tpose = 0.0_rk
    no2_xs_tpose= 0.0_rk

    call o3xs( nlev,tlev,nwave,wl,o3_xs_tpose )
    call no2xs_jpl06a( nlev,tlev,nwave,wl,no2_xs_tpose )
    o3_xs  = transpose( o3_xs_tpose )
    no2_xs = transpose( no2_xs_tpose )

    call tuv_radfld( nlambda_start, cld_od_opt, cldfrac, nlyr, nwave, &
         zen, zlev, alb, &
         aircol, o3col, so2col, no2col, &
         tauaer300, tauaer400, tauaer600, tauaer999, &
         waer300, waer400, waer600, waer999, &
         gaer300, gaer400, gaer600, gaer999, &
         dtaer, omaer, gaer, dtcld, omcld, gcld, &
         has_aer_ra_feedback, &
         qll, dobsi, o3_xs, no2_xs, o2_xs, &
         so2_xs, wl(1), wc, tlev, dto2, radfld, efld, &
         e_dir, e_dn, e_up, &
         dir_fld, dwn_fld, up_fld, dt_cld, errmsg, errflg )

  end subroutine tuv_radiation_transfer_run
  
!> \section arg_table_tuv_radiation_transfer_finalize Argument Table
!! | local_name | standard_name                         | long_name                      | units     | rank | type      | kind      | intent | optional |
!! |------------|---------------------------------------|--------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message                    | CCPP error message             | none      |    0 | character | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag                       | CCPP error flag                | flag      |    0 | integer   |           | out    | F        |
!!
  subroutine tuv_radiation_transfer_finalize( errmsg, errflg )

    !--- arguments
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine tuv_radiation_transfer_finalize

end module tuv_radiation_transfer
