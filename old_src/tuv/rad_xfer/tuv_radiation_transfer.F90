module tuv_radiation_transfer
  use phot_kind_mod, only: rk => kind_phot
  use phot_kind_mod, only: kind_phys => kind_phot

  implicit none

  private
  public :: tuv_radiation_transfer_init
  public :: tuv_radiation_transfer_run
  
contains
  
!> \section arg_table_tuv_radiation_transfer_init Argument Table
!! \htmlinclude tuv_radiation_transfer_init.html
!!
subroutine tuv_radiation_transfer_init( realkind, errmsg, errflg )
    use wavelength_grid,  only: nwave, wl
    use module_xsections, only: rdxs_init

    integer,          intent(in)  :: realkind
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    
    errmsg = ''
    errflg = 0

    if ( realkind/=kind_phys ) then
       errmsg = 'tuv_radiation_transfer_init: realkind does not match kind_phot'
       errflg = 1
       return
    end if
    
    call rdxs_init( nwave, wl, errmsg, errflg )
    
    if (errflg.ne.0) return
   
  end subroutine tuv_radiation_transfer_init
  
!> \section arg_table_tuv_radiation_transfer_run Argument Table
!! \htmlinclude tuv_radiation_transfer_run.html
!!
  subroutine tuv_radiation_transfer_run( nlayer, tuv_n_wavelen, zenith, albedo, press_mid, press_top, alt, temp, &
       o3vmr, so2vmr, no2vmr, cldfrc, cldwat, dto2, radfld, errmsg, errflg )

    use tuv_subs,         only: tuv_radfld
    use wavelength_grid,  only: nwave, wl, wc
    use module_xsections, only: o2_xs, so2_xs, o3xs, no2xs_jpl06a
    use params_mod
 
    integer,          intent(in)  :: nlayer

    !! NOTE THIS VARIABLE WILL GO AWAY - FOR NOW IS REQUIRED WORKAROUND FOR CPF
    integer,          intent(in)    :: tuv_n_wavelen

    real(kind_phys),  intent(in)  :: zenith
    real(kind_phys),  intent(in)  :: albedo
    real(kind_phys),  intent(in)  :: press_mid(:)
    real(kind_phys),  intent(in)  :: press_top
    real(kind_phys),  intent(in)  :: alt(:)  ! m
    real(kind_phys),  intent(in)  :: temp(:) ! K
    real(kind_phys),  intent(in)  :: o3vmr(:)
    real(kind_phys),  intent(in)  :: so2vmr(:)
    real(kind_phys),  intent(in)  :: no2vmr(:)
    real(kind_phys),  intent(in)  :: cldfrc(:)
    real(kind_phys),  intent(in)  :: cldwat(:) ! cld water content (g/m3)
    real(kind_phys),  intent(in)  :: dto2(:,:)
    real(kind_phys),  intent(out) :: radfld(:,:) ! /sec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: nlambda_start=1
    integer, parameter :: cld_od_opt=1
    logical, parameter :: has_aer_ra_feedback = .false.
    real(rk), parameter :: dobsi = 0._rk
    
    real(rk) :: zen
    real(rk) :: alb(nwave)
    real(rk) :: zlev(nlayer+1) ! km 
    real(rk) :: tlev(nlayer)
    real(rk) :: cldfrclev(nlayer)
    real(rk) :: cldwatlev(nlayer)
    real(rk) :: aircol(nlayer)  ! # molecules / cm2 in each layer
    real(rk) :: o3col(nlayer) 
    real(rk) :: so2col(nlayer)
    real(rk) :: no2col(nlayer)
    real(rk) :: dpress(nlayer)

    real(rk) :: tauaer300(nlayer) ! aerosol properties
    real(rk) :: tauaer400(nlayer)
    real(rk) :: tauaer600(nlayer)
    real(rk) :: tauaer999(nlayer)
    real(rk) :: waer300(nlayer)
    real(rk) :: waer400(nlayer)
    real(rk) :: waer600(nlayer)
    real(rk) :: waer999(nlayer)
    real(rk) :: gaer300(nlayer)
    real(rk) :: gaer400(nlayer)
    real(rk) :: gaer600(nlayer)
    real(rk) :: gaer999(nlayer)
    
    real(rk) :: dtaer(nlayer,nwave), omaer(nlayer,nwave), gaer(nlayer,nwave)
    real(rk) :: dtcld(nlayer,nwave), omcld(nlayer,nwave), gcld(nlayer,nwave)
    real(rk) :: dt_cld(nlayer)
    
    real(rk) :: efld(nlayer,nwave)
    real(rk) :: e_dir(nlayer,nwave)
    real(rk) :: e_dn(nlayer,nwave)
    real(rk) :: e_up(nlayer,nwave)
    real(rk) :: dir_fld(nlayer,nwave)
    real(rk) :: dwn_fld(nlayer,nwave)
    real(rk) :: up_fld(nlayer,nwave)

    integer :: k, kk

    real(rk) :: o3_xs(nwave,nlayer)
    real(rk) :: no2_xs(nwave,nlayer)
    real(rk) :: o3_xs_tpose(nlayer,nwave)
    real(rk) :: no2_xs_tpose(nlayer,nwave)
    real(rk) :: delz_km, delz_cm
    
    errmsg = ''
    errflg = 0

    dpress(nlayer-1:1:-1) = press_mid(2:nlayer) - press_mid(1:nlayer-1)
    do k=1,nlayer-1
       kk=nlayer-1-k+1
       aircol(k) = 10._rk*dpress(k)*R/(kboltz*g)
       o3col(k)  = 0.5_rk*(o3vmr(kk)+o3vmr(kk+1))*aircol(k)
       so2col(k) = 0.5_rk*(so2vmr(kk)+so2vmr(kk+1))*aircol(k)
       no2col(k) = 0.5_rk*(no2vmr(kk)+no2vmr(kk+1))*aircol(k)
    end do

    ! inputs need to be bottom up vert coord
    tlev(nlayer:1:-1) = temp(1:nlayer)
    cldfrclev(nlayer:1:-1) = cldfrc(1:nlayer)
    cldwatlev(nlayer:1:-1) = cldwat(1:nlayer)
    zlev(nlayer:1:-1) = alt(1:nlayer)*1.e-3_rk ! m -> km

    delz_km = zlev(nlayer) - zlev(nlayer-1)
    delz_cm = delz_km*1.e5_rk ! cm

    zlev(nlayer+1) = zlev(nlayer) + delz_km ! km
    aircol(nlayer) = delz_cm * 10._rk * press_top / ( kboltz * tlev(nlayer) ) ! molecules / cm2
    o3col(nlayer)  = o3vmr(1)  * aircol(nlayer)
    so2col(nlayer) = so2vmr(1) * aircol(nlayer)
    no2col(nlayer) = no2vmr(1) * aircol(nlayer)

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

    call o3xs( nlayer,tlev,nwave,wl,o3_xs_tpose )
    call no2xs_jpl06a( nlayer,tlev,nwave,wl,no2_xs_tpose )
    o3_xs  = transpose( o3_xs_tpose )
    no2_xs = transpose( no2_xs_tpose )

    call tuv_radfld( nlambda_start, cld_od_opt, cldfrclev, nlayer, nwave, &
         zen, zlev, alb, &
         aircol, o3col, so2col, no2col, &
         tauaer300, tauaer400, tauaer600, tauaer999, &
         waer300, waer400, waer600, waer999, &
         gaer300, gaer400, gaer600, gaer999, &
         dtaer, omaer, gaer, dtcld, omcld, gcld, &
         has_aer_ra_feedback, &
         cldwatlev, dobsi, o3_xs, no2_xs, o2_xs, &
         so2_xs, wl(1), wc, tlev, dto2, radfld, efld, &
         e_dir, e_dn, e_up, &
         dir_fld, dwn_fld, up_fld, dt_cld, errmsg, errflg )

  end subroutine tuv_radiation_transfer_run

end module tuv_radiation_transfer
