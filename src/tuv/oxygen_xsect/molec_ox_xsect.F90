module molec_ox_xsect
  use phot_kind_mod, only: rk => kind_phot
  use phot_kind_mod, only: kind_phys => kind_phot
  use wavelength_grid, only: nwave, wl

  implicit none

contains

!> \section arg_table_molec_ox_xsect_init Argument Table
!! \htmlinclude molec_ox_xsect_init.html
!!
  subroutine molec_ox_xsect_init( errmsg, errflg )
    use la_srb_mod, only : la_srb_init

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call la_srb_init( errmsg, errflg )

  end subroutine molec_ox_xsect_init

!> \section arg_table_molec_ox_xsect_run Argument Table
!! \htmlinclude molec_ox_xsect_run.html
!!
  subroutine molec_ox_xsect_run( nlev, zen, alt, temp, press_mid, press_top, o2vmr, dto2, srb_o2_xs, errmsg, errflg )
    use module_xsections, only: o2_xs
    use phot_util_mod, only : sphers, airmas
    use la_srb_mod,    only : la_srb_comp
    use params_mod,    only : R, g, kboltz

    integer,          intent(in)    :: nlev
    real(kind_phys),  intent(in)    :: zen
    real(kind_phys),  intent(in)    :: alt(:)  ! m
    real(kind_phys),  intent(in)    :: temp(:) ! K
    real(kind_phys),  intent(in)    :: press_mid(:)
    real(kind_phys),  intent(in)    :: press_top
    real(kind_phys),  intent(in)    :: o2vmr(:)
    real(kind_phys),  intent(out)   :: dto2(:,:)
    real(kind_phys),  intent(out)   :: srb_o2_xs(:,:)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer  :: nlyr, k
    integer  :: wn
    integer  :: nid(0:nlev)
    real(rk) :: dsdh(0:nlev,nlev)
    real(rk) :: vcol(nlev)
    real(rk) :: scol(nlev)
    real(rk) :: tlev(nlev)
    real(rk) :: zlev(nlev+1)
    real(rk) :: o2lev(nlev+1)

    real(rk) :: aircol(nlev)  ! # molecules / cm2 in each layer
    real(rk) :: o2col(nlev)  
    real(rk) :: dpress(nlev-1)
    real(rk) :: delz_cm, delz_km

    errmsg = ''
    errflg = 0

    dto2(:,:) = 0._rk
    srb_o2_xs(:,:) = 0._rk

    nlyr = nlev-1

    dpress(nlyr:1:-1) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    zlev(nlev:1:-1) = alt(1:nlev) *1.e-3_rk ! m -> km
    o2lev(nlev:1:-1) = o2vmr(1:nlev)

    delz_km = zlev(nlev) - zlev(nlev-1)
    delz_cm = delz_km*1.e5_rk ! cm

    zlev(nlev+1) = zlev(nlev) + delz_km ! km

    do k=1,nlyr
       aircol(k) = 10._rk*dpress(k)*R/(kboltz*g)
       o2col(k)  = 0.5_rk*(o2lev(k)+o2lev(k+1))*aircol(k)
    end do

    aircol(nlev) = delz_cm * 10._rk * press_top / ( kboltz * temp(1) ) ! molecules / cm2
    o2col(nlev)  = o2lev(nlev) * aircol(nlev)

    tlev(nlev:1:-1) = temp(1:nlev)

    do wn = 1,nwave
       dto2(:nlev,wn) = o2col(:nlev) * o2_xs(wn)
    end do

    call sphers( nlev, zlev, zen, dsdh, nid )
    call airmas( nlev, dsdh, nid, aircol, vcol, scol )

    call la_srb_comp( nlev, wl(1), tlev, vcol, scol, o2lev, o2_xs, dto2, srb_o2_xs, errmsg, errflg )

  end subroutine molec_ox_xsect_run

end module molec_ox_xsect
