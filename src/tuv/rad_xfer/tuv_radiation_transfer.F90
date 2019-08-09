module tuv_radiation_transfer
!  use phot_kind_mod, only: rk => kind_phot
  use phot_kind_mod, only: kind_phys => kind_phot

  implicit none

  private
  public :: tuv_radiation_transfer_init
  public :: tuv_radiation_transfer_run
  public :: tuv_radiation_transfer_finalize
  
  integer :: nlev, nlyr
  
contains
  
!> \section arg_table_tuv_radiation_transfer_init Argument Table
!! \htmlinclude tuv_radiation_transfer_init.html
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

    if ( realkind/=kind_phys ) then
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
!! \htmlinclude tuv_radiation_transfer_run.html
!!
  subroutine tuv_radiation_transfer_run( nlev, tuv_n_wavelen, zenith, albedo, press_mid, alt, temp, o3vmr, &
         so2vmr, no2vmr, dto2, radfld, errmsg, errflg )

    use tuv_subs,         only: tuv_radfld
    use rad_abs_xsect,    only: o2_xs, so2_xs, nwave, wl, wc
    use module_xsections, only: o3xs, no2xs_jpl06a
    use params_mod
 
    integer,          intent(in)  :: nlev
    integer,          intent(in)  :: tuv_n_wavelen
    real(kind_phys),         intent(in)  :: zenith
    real(kind_phys),         intent(in)  :: albedo
    real(kind_phys),         intent(in)  :: press_mid(:)
    real(kind_phys),         intent(in)  :: alt(:)  ! m
    real(kind_phys),         intent(in)  :: temp(:) ! K
    real(kind_phys),         intent(in)  :: o3vmr(:)
    real(kind_phys),         intent(in)  :: so2vmr(:)
    real(kind_phys),         intent(in)  :: no2vmr(:)
    real(kind_phys),         intent(in)  :: dto2(:,:)
    real(kind_phys),         intent(out) :: radfld(:,:) ! /sec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: nlambda_start=1
    integer, parameter :: cld_od_opt=1
    logical, parameter :: has_aer_ra_feedback = .false.
    real(kind_phys), parameter :: dobsi = 0._kind_phys
    
    real(kind_phys) :: zen
    real(kind_phys) :: alb(nwave)
    real(kind_phys) :: zlev(nlev) ! km 
    real(kind_phys) :: tlev(nlev)
    real(kind_phys) :: aircol(nlyr)  ! # molecules / cm2 in each layer
    real(kind_phys) :: o3col(nlyr) 
    real(kind_phys) :: so2col(nlyr)
    real(kind_phys) :: no2col(nlyr)
    real(kind_phys) :: dpress(nlyr)

    real(kind_phys) :: tauaer300(nlev) ! aerosol properties
    real(kind_phys) :: tauaer400(nlev)
    real(kind_phys) :: tauaer600(nlev)
    real(kind_phys) :: tauaer999(nlev)
    real(kind_phys) :: waer300(nlev)
    real(kind_phys) :: waer400(nlev)
    real(kind_phys) :: waer600(nlev)
    real(kind_phys) :: waer999(nlev)
    real(kind_phys) :: gaer300(nlev)
    real(kind_phys) :: gaer400(nlev)
    real(kind_phys) :: gaer600(nlev)
    real(kind_phys) :: gaer999(nlev)
    
    real(kind_phys) :: dtaer(nlyr,nwave), omaer(nlyr,nwave), gaer(nlyr,nwave)
    real(kind_phys) :: dtcld(nlyr,nwave), omcld(nlyr,nwave), gcld(nlyr,nwave)
    real(kind_phys) :: dt_cld(nlyr)
    
    real(kind_phys) :: qll(nlev) ! cld water content (g/m3)
    real(kind_phys) :: cldfrac(nlev)
    real(kind_phys) :: efld(nlev,nwave)
    real(kind_phys) :: e_dir(nlev,nwave)
    real(kind_phys) :: e_dn(nlev,nwave)
    real(kind_phys) :: e_up(nlev,nwave)
    real(kind_phys) :: dir_fld(nlev,nwave)
    real(kind_phys) :: dwn_fld(nlev,nwave)
    real(kind_phys) :: up_fld(nlev,nwave)

    integer :: k, kk

    real(kind_phys) :: o3_xs(nwave,nlev)
    real(kind_phys) :: no2_xs(nwave,nlev)
    real(kind_phys) :: o3_xs_tpose(nlev,nwave)
    real(kind_phys) :: no2_xs_tpose(nlev,nwave)

    errmsg = ''
    errflg = 0

    dpress(1:nlyr) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    do k=1,nlyr
       kk=nlyr-k+1
       aircol(k) = 10._kind_phys*dpress(k)*R/(kboltz*g)
       o3col(kk)  = 0.5_kind_phys*(o3vmr(k)+o3vmr(k+1))*aircol(k)
       so2col(kk) = 0.5_kind_phys*(so2vmr(k)+so2vmr(k+1))*aircol(k)
       no2col(kk) = 0.5_kind_phys*(no2vmr(k)+no2vmr(k+1))*aircol(k)
    end do

    ! inputs need to be bottom up vert coord
    aircol(1:nlyr) = aircol(nlyr:1:-1)
    tlev(nlev:1:-1) = temp(1:nlev)
    zlev(nlev:1:-1) = alt(1:nlev)*1.e-3_kind_phys ! m -> km

    qll=0.0_kind_phys
    cldfrac=0.0_kind_phys
    tauaer300=0.0_kind_phys
    tauaer400=0.0_kind_phys
    tauaer600=0.0_kind_phys
    tauaer999=0.0_kind_phys
    waer300=1.0_kind_phys
    waer400=1.0_kind_phys
    waer600=1.0_kind_phys
    waer999=1.0_kind_phys
    gaer300=0.0_kind_phys
    gaer400=0.0_kind_phys
    gaer600=0.0_kind_phys
    gaer999=0.0_kind_phys

    zen = zenith
    alb(:) = albedo

    o3_xs_tpose = 0.0_kind_phys
    no2_xs_tpose= 0.0_kind_phys

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
!! \htmlinclude tuv_radiation_transfer_finalize.html
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
