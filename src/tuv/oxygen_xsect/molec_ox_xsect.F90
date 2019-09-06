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
  subroutine molec_ox_xsect_run( nlyr, zen, alt, temp, press_mid, press_top, o2vmr, dto2, srb_o2_xs, errmsg, errflg )
    use module_xsections, only: o2_xs
    use phot_util_mod, only : sphers, airmas
    use la_srb_mod,    only : la_srb_comp
    use params_mod,    only : R, g, kboltz

    integer,          intent(in)    :: nlyr
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

    integer  :: k
    integer  :: wn
    integer  :: nid(0:nlyr)
    real(rk) :: dsdh(0:nlyr,nlyr)
    real(rk) :: vcol(nlyr)
    real(rk) :: scol(nlyr)
    real(rk) :: tlev(nlyr)
    real(rk) :: zlev(nlyr+1)
    real(rk) :: o2lev(nlyr+1)

    real(rk) :: aircol(nlyr)  ! # molecules / cm2 in each layer
    real(rk) :: o2col(nlyr)  
    real(rk) :: dpress(nlyr-1)
    real(rk) :: delz_cm, delz_km

    errmsg = ''
    errflg = 0

    dto2(:,:) = 0._rk
    srb_o2_xs(:,:) = 0._rk

    dpress(nlyr-1:1:-1) = press_mid(2:nlyr) - press_mid(1:nlyr-1)
    zlev(nlyr:1:-1) = alt(1:nlyr) *1.e-3_rk ! m -> km
    o2lev(nlyr:1:-1) = o2vmr(1:nlyr)

    delz_km = zlev(nlyr) - zlev(nlyr-1)
    delz_cm = delz_km*1.e5_rk ! cm

    zlev(nlyr+1) = zlev(nlyr) + delz_km ! km

    do k=1,nlyr-1
       aircol(k) = 10._rk*dpress(k)*R/(kboltz*g)
       o2col(k)  = 0.5_rk*(o2lev(k)+o2lev(k+1))*aircol(k)
    end do

    aircol(nlyr) = delz_cm * 10._rk * press_top / ( kboltz * temp(1) ) ! molecules / cm2
    o2col(nlyr)  = o2lev(nlyr) * aircol(nlyr)

    tlev(nlyr:1:-1) = temp(1:nlyr)

    do wn = 1,nwave
       dto2(:nlyr,wn) = o2col(:nlyr) * o2_xs(wn)
    end do

    call sphers( nlyr, zlev, zen, dsdh, nid )
    call airmas( nlyr, dsdh, nid, aircol, vcol, scol )

    call la_srb_comp( nlyr, wl(1), tlev, vcol, scol, o2lev, o2_xs, dto2, srb_o2_xs, errmsg, errflg )

  end subroutine molec_ox_xsect_run

end module molec_ox_xsect
