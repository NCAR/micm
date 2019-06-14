module molec_ox_xsect
  use phot_kind_mod, only: rk => kind_phot
  use rad_abs_xsect, only: o2_xs, nwave, wl

  implicit none

contains

!> \section arg_table_molec_ox_xsect_init Argument Table
!! | local_name | standard_name             | long_name                 | units   | rank | type      | kind      | intent | optional |
!! |------------|---------------------------|---------------------------|---------|------|-----------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message        | CCPP error message        | none    |    0 | character | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag           | CCPP error flag           | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine molec_ox_xsect_init( errmsg, errflg )
    use la_srb_mod, only : la_srb_init

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call la_srb_init( errmsg, errflg )

  end subroutine molec_ox_xsect_init

!> \section arg_table_molec_ox_xsect_run Argument Table
!! | local_name | standard_name             | long_name                          | units     | rank | type        | kind      | intent | optional |
!! |------------|---------------------------|------------------------------------|-----------|------|-------------|-----------|--------|----------|
!! | nlev       | num_levels_for_photolysis | number of column layers            | count     |    0 | integer     |           | in     | F        |
!! | zen        | solar_zenith              | solar zenith angle                 | degrees   |    0 | real        | kind_phys | in     | F        |
!! | alt        | layer_altitude            | mid-point layer altitude           | km        |    1 | real        | kind_phys | in     | F        |
!! | temp       | layer_temperature         | mid-point layer temperature        | K         |    1 | real        | kind_phys | in     | F        |
!! | press_mid  | layer_pressure            | mid-point layer pressure           | Pa        |    1 | real        | kind_phys | in     | F        |
!! | o2vmr      | O2_vmr_col                | O2 volume mixing ratio column      | mole/mole |    1 | real        | kind_phys | in     | F        |
!! | dto2       | O2_optical_depth          | optical depth due to O2 absorption | cm        |    2 | real        | kind_phys | in     | F        |
!! | srb_o2_xs  | O2_xsect                  | O2 effective cross section         | cm2       |    2 | real        | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message        | CCPP error message                 | none      |    0 | character   | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag           | CCPP error flag                    | flag      |    0 | integer     |           | out    | F        |
!!
  subroutine molec_ox_xsect_run( nlev, zen, alt, temp, press_mid, o2vmr, dto2, srb_o2_xs, errmsg, errflg )
    use rad_abs_xsect, only : o2_xs
    use phot_util_mod, only : sphers, airmas
    use la_srb_mod,    only : la_srb_comp
    use params_mod,    only : R, g, kboltz

    integer,          intent(in)    :: nlev
    real(rk),         intent(in)    :: zen
    real(rk),         intent(in)    :: alt(:)  ! km
    real(rk),         intent(in)    :: temp(:) ! K
    real(rk),         intent(in)    :: press_mid(:)
    real(rk),         intent(in)    :: o2vmr(:)
    real(rk),         intent(inout) :: dto2(:,:)
    real(rk),         intent(inout) :: srb_o2_xs(:,:)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer  :: nlyr, k
    integer  :: wn
    integer  :: nid(0:nlev-1)
    real(rk) :: dsdh(0:nlev-1,nlev-1)
    real(rk) :: vcol(nlev-1)
    real(rk) :: scol(nlev-1)
    real(rk) :: tlev(nlev)
    real(rk) :: zlev(nlev)

    real(rk) :: aircol(nlev-1)  ! # molecules / cm2 in each layer
    real(rk) :: o2col(nlev)  
    real(rk) :: dpress(nlev)

    errmsg = ''
    errflg = 0

    dto2(:,:) = 0._rk
    srb_o2_xs(:,:) = 0._rk

    nlyr = nlev-1
    dpress(1:nlyr) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    do k=1,nlyr
       aircol(k) = 10._rk*dpress(k)*R/(kboltz*g)
       o2col(k)  = 0.5_rk*(o2vmr(k)+o2vmr(k+1))*aircol(k)
    end do
    aircol(1:nlyr) = aircol(nlyr:1:-1)
    o2col(1:nlyr)  = o2col(nlyr:1:-1)
    tlev(nlev:1:-1) = temp(1:nlev)
    zlev(nlev:1:-1) = alt(1:nlev) ! km

    do wn = 1,nwave
       dto2(:nlyr,wn) = o2col(:nlyr) * o2_xs(wn)
    end do

    call sphers( nlyr, zlev, zen, dsdh, nid )
    call airmas( nlyr, dsdh, nid, aircol, vcol, scol )

    call la_srb_comp( nlyr, wl(1), tlev, vcol, scol, o2vmr, o2_xs, dto2, srb_o2_xs )

  end subroutine molec_ox_xsect_run

!> \section arg_table_molec_ox_xsect_finalize Argument Table
!! | local_name | standard_name      | long_name          | units     | rank | type      | kind  | intent | optional |
!! |------------|--------------------|--------------------|-----------|------|-----------|-------|--------|----------|
!! | errmsg     | ccpp_error_message | CCPP error message | none      |    0 | character | len=* | out    | F        |
!! | errflg     | ccpp_error_flag    | CCPP error flag    | flag      |    0 | integer   |       | out    | F        |
!!
  subroutine molec_ox_xsect_finalize( errmsg, errflg )

    !--- arguments
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine molec_ox_xsect_finalize

end module molec_ox_xsect
