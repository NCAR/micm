module tuv_photolysis

  use phot_kind_mod, only: rk => kind_phot
  use module_prates_tuv, only: calc_tuv_init, calc_tuv_prates
  use params_mod, only: input_data_root

  implicit none

  integer, protected :: tuv_n_wavelen = -1
  integer, parameter :: tuv_n_phot = 113
  character(len=16), parameter :: tuv_jnames(tuv_n_phot) = &
       (/'j_o2            ' &
       , 'j_o1d           ' &
       , 'j_o3p           ' &
       , 'j_no2           ' &
       , 'j_no3_a         ' &
       , 'j_no3_b         ' &
       , 'j_n2o5_a        ' &
       , 'j_n2o5_b        ' &
       , 'j_hno2          ' &
       , 'j_hno3          ' &
       , 'j_hno4          ' &
       , 'j_h2o2          ' &
       , 'j_chbr3         ' &
       , 'j_ch3cho_a      ' &
       , 'j_ch3cho_b      ' &
       , 'j_ch3cho_c      ' &
       , 'j_c2h5cho       ' &
       , 'j_gly_a         ' &
       , 'j_gly_b         ' &
       , 'j_gly_c         ' &
       , 'j_mgly          ' &
       , 'j_ch3coch3      ' &
       , 'j_ch3ooh        ' &
       , 'j_ch3ono2       ' &
       , 'j_pan_a         ' &
       , 'j_pan_b         ' &
       , 'j_ccl2o         ' &
       , 'j_ccl4          ' &
       , 'j_cclfo         ' &
       , 'j_cf2o          ' &
       , 'j_cf2clcfcl2    ' &
       , 'j_cf2clcf2cl    ' &
       , 'j_cf3cf2cl      ' &
       , 'j_ccl3f         ' &
       , 'j_ccl2f2        ' &
       , 'j_ch3br         ' &
       , 'j_ch3ccl3       ' &
       , 'j_ch3cl         ' &
       , 'j_cloo          ' &
       , 'j_cf3chcl2      ' &
       , 'j_cf3chfcl      ' &
       , 'j_ch3cfcl2      ' &
       , 'j_ch3cf2cl      ' &
       , 'j_cf3cf2chcl2   ' &
       , 'j_cf2clcf2chfcl ' &
       , 'j_chclf2        ' &
       , 'j_ho2           ' &
       , 'j_cf2bf2        ' &
       , 'j_cf2brcl       ' &
       , 'j_cf3br         ' &
       , 'j_cf2brcf2br    ' &
       , 'j_n2o           ' &
       , 'j_clono2_a      ' &
       , 'j_clono2_b      ' &
       , 'j_brono2_a      ' &
       , 'j_brono2_b      ' &
       , 'j_cl2           ' &
       , 'j_glyald_a      ' &
       , 'j_glyald_b      ' &
       , 'j_glyald_c      ' &
       , 'j_biacetyl      ' &
       , 'j_mvk           ' &
       , 'j_macr          ' &
       , 'j_ch3cocooh     ' &
       , 'j_ch3ch2ono2    ' &
       , 'j_ch3chono2ch3  ' &
       , 'j_ch2ohch2ono2  ' &
       , 'j_ch3coch2ono2  ' &
       , 'j_bnit1         ' &
       , 'j_cloocl        ' &
       , 'j_hyac_a        ' &
       , 'j_hyac_b        ' &
       , 'j_hobr          ' &
       , 'j_bro           ' &
       , 'j_br2           ' &
       , 'j_no3_aq_a      ' &
       , 'j_no3_aq_b      ' &
       , 'j_no3_aq_c      ' &
       , 'j_mek           ' &
       , 'j_ppn_a         ' &
       , 'j_ppn_b         ' &
       , 'j_hoch2ooh      ' &
       , 'j_acrol         ' &
       , 'j_ch3coooh      ' &
       , 'j_amine         ' &
       , 'j_clo_a         ' &
       , 'j_clo_b         ' &
       , 'j_clno2         ' &
       , 'j_brno          ' &
       , 'j_brno2         ' &
       , 'j_brono_a       ' &
       , 'j_brono_b       ' &
       , 'j_hocl          ' &
       , 'j_nocl          ' &
       , 'j_oclo          ' &
       , 'j_brcl          ' &
       , 'j_ch3oono2      ' &
       , 'j_bnit2         ' &
       , 'j_clono         ' &
       , 'j_hcl           ' &
       , 'j_ch2o_r        ' &
       , 'j_ch2o_m        ' &
       , 'j_ch3cooh       ' &
       , 'j_ch3ocl        ' &
       , 'j_chcl3         ' &
       , 'j_c2h5ono2      ' &
       , 'j_nc3h7ono2     ' &
       , 'j_1c4h9ono2     ' &
       , 'j_2c4h9ono2     ' &
       , 'j_perfluoro     ' &
       , 'j_i2            ' &
       , 'j_io            ' &
       , 'j_ioh           ' &
       /)
  
contains

  subroutine tuv_photolysis_readnl(nml_file) ! this will be a CPF interface someday

    use module_prates_tuv, only: get_xsqy_tab, nwave
 
    character(len=*), intent(in)  :: nml_file

    character(len=512) :: errmsg
    integer :: errflg

    character(len=512) :: xsqy_filepath

    namelist /tuv_opts/ input_data_root

    open(unit=10,file=nml_file)
    read(unit=10,nml=tuv_opts)
    close(10)

    xsqy_filepath = trim(input_data_root)//'/wrf_tuv_xsqy.nc'
    call get_xsqy_tab(xsqy_filepath, errmsg, errflg) ! call this here since nwave needs to be known earlier than the init phase
    tuv_n_wavelen = nwave
    
  end subroutine tuv_photolysis_readnl

!> \section arg_table_tuv_photolysis_init Argument Table
!! | local_name | standard_name             | long_name                 | units   | rank | type      | kind      | intent | optional |
!! |------------|---------------------------|---------------------------|---------|------|-----------|-----------|--------|----------|
!! | realkind   | phys_real_kind            | physics real kind         | none    |    0 | integer   |           | in     | F        |
!! | errmsg     | ccpp_error_message        | CCPP error message        | none    |    0 | character | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag           | CCPP error flag           | flag    |    0 | integer   |           | out    | F        |
!!
subroutine tuv_photolysis_init( realkind, errmsg, errflg )

    integer,          intent(in)  :: realkind
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    logical, parameter :: full_tuv = .true.

    errmsg = ' '
    errflg = 0

    if ( realkind/=rk ) then
       errmsg = 'tuv_photolysis_init: realkind does not match kind_phot'
       errflg = 1
       return
    end if

    call  calc_tuv_init( full_tuv, tuv_jnames, errmsg, errflg )

  end subroutine tuv_photolysis_init

!> \section arg_table_tuv_photolysis_run Argument Table
!! | local_name | standard_name                         | long_name                      | units     | rank | type      | kind      | intent | optional |
!! |------------|---------------------------------------|--------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | nlev       | num_levels_for_photolysis             | number of column layers        | count     |    0 | integer   |           | in     | F        |
!! | temp       | layer_temperature                     | mid-point layer temperature    | K         |    1 | real      | kind_phys | in     | F        |
!! | press_mid  | layer_pressure                        | mid-point layer pressure       | Pa        |    1 | real      | kind_phys | in     | F        |
!! | radfld     | actinic_photon_fluxes                 | actinic photon fluxes          | cm-2 sec-1|    2 | real      | kind_phys | in     | F        |
!! | srb_o2_xs  | O2_xsect                              | O2 cross sections              | cm2       |    2 | real      | kind_phys | in     | F        |
!! | tuv_prates | photolysis_rates_col                  | photolysis rates column        | s-1       |    2 | real      | kind_phys | out    | F        |
!! | errmsg     | ccpp_error_message                    | CCPP error message             | none      |    0 | character | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag                       | CCPP error flag                | flag      |    0 | integer   |           | out    | F        |
!!
subroutine tuv_photolysis_run( nlev, temp, press_mid, radfld, srb_o2_xs, tuv_prates, errmsg, errflg )

    integer,          intent(in)  :: nlev
    real(rk),         intent(in)  :: temp(:)
    real(rk),         intent(in)  :: press_mid(:)
    real(rk),         intent(in)  :: radfld(:,:) ! (nwave,nlev)
    real(rk),         intent(in)  :: srb_o2_xs(:,:) !(nwave,kts:kte)
    real(rk),         intent(out) :: tuv_prates(:,:) ! /sec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k, kk, j
    real(rk) :: airdens(nlev) ! # molecules / cm3 in each layer
    real(rk) :: tlev(nlev) ! # K -- bottom up

    real(rk), parameter :: kboltz= 1.38064852e-16_rk ! boltzmann constant (erg/K)

    ! inputs need to be bottom vertical coord
    do k=1,nlev
       kk=nlev-k+1
       airdens(kk) = 10._rk*press_mid(k)/(kboltz*temp(k))
    end do
    tlev(nlev:1:-1) = temp(1:nlev)

    call calc_tuv_prates(1 ,nlev,nlev, tlev, airdens, radfld, srb_o2_xs, tuv_prates, errmsg, errflg)

    ! return top down rates
    do j=1,tuv_n_phot
       tuv_prates(:nlev,j) = tuv_prates(nlev:1:-1,j)
    end do
    
  end subroutine tuv_photolysis_run
  
!> \section arg_table_tuv_photolysis_finalize Argument Table
!! | local_name | standard_name                         | long_name                      | units     | rank | type      | kind      | intent | optional |
!! |------------|---------------------------------------|--------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message                    | CCPP error message             | none      |    0 | character | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag                       | CCPP error flag                | flag      |    0 | integer   |           | out    | F        |
!!
  subroutine tuv_photolysis_finalize( errmsg, errflg )

    !--- arguments
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine tuv_photolysis_finalize

end module tuv_photolysis
