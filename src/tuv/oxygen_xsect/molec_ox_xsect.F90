module molec_ox_xsect
!  use phot_kind_mod, only: rk => kind_phot
  use phot_kind_mod, only: kind_phys => kind_phot
  use rad_abs_xsect, only: o2_xs, nwave, wl

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
  subroutine molec_ox_xsect_run( nlev, tuv_n_wavelen, nlevelsMinus1, zen, alt, temp, press_mid, o2vmr, dto2, srb_o2_xs, &
         errmsg, errflg )
    use rad_abs_xsect, only : o2_xs
    use phot_util_mod, only : sphers, airmas
    use la_srb_mod,    only : la_srb_comp
    use params_mod,    only : R, g, kboltz

    integer,          intent(in)    :: nlev
    integer,          intent(in)    :: tuv_n_wavelen

    !! NOTE THIS VARIABLE WILL GO AWAY - FOR NOW IS REQUIRED WORKAROUND FOR CPF
    integer,          intent(in)    :: nlevelsMinus1

    real(kind_phys),         intent(in)    :: zen
    real(kind_phys),         intent(in)    :: alt(:)  ! m
    real(kind_phys),         intent(in)    :: temp(:) ! K
    real(kind_phys),         intent(in)    :: press_mid(:)
    real(kind_phys),         intent(in)    :: o2vmr(:)
    real(kind_phys),         intent(out) :: dto2(:,:)
    real(kind_phys),         intent(out) :: srb_o2_xs(:,:)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer  :: nlyr, k
    integer  :: wn
    integer  :: nid(0:nlev-1)
    real(kind_phys) :: dsdh(0:nlev-1,nlev-1)
    real(kind_phys) :: vcol(nlev-1)
    real(kind_phys) :: scol(nlev-1)
    real(kind_phys) :: tlev(nlev)
    real(kind_phys) :: zlev(nlev)

    real(kind_phys) :: aircol(nlev-1)  ! # molecules / cm2 in each layer
    real(kind_phys) :: o2col(nlev)  
    real(kind_phys) :: dpress(nlev)

    errmsg = ''
    errflg = 0

    dto2(:,:) = 0._kind_phys
    srb_o2_xs(:,:) = 0._kind_phys

    nlyr = nlev-1
    dpress(1:nlyr) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    do k=1,nlyr
       aircol(k) = 10._kind_phys*dpress(k)*R/(kboltz*g)
       o2col(k)  = 0.5_kind_phys*(o2vmr(k)+o2vmr(k+1))*aircol(k)
    end do
    aircol(1:nlyr) = aircol(nlyr:1:-1)
    o2col(1:nlyr)  = o2col(nlyr:1:-1)
    tlev(nlev:1:-1) = temp(1:nlev)
    zlev(nlev:1:-1) = alt(1:nlev)*1.e-3_kind_phys ! m -> km

    do wn = 1,nwave
       dto2(:nlyr,wn) = o2col(:nlyr) * o2_xs(wn)
    end do

    call sphers( nlyr, zlev, zen, dsdh, nid )
    call airmas( nlyr, dsdh, nid, aircol, vcol, scol )

    call la_srb_comp( nlyr, wl(1), tlev, vcol, scol, o2vmr, o2_xs, dto2, srb_o2_xs )

  end subroutine molec_ox_xsect_run

end module molec_ox_xsect
