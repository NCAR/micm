module tuv_photolysis

  use phot_kind_mod, only: kind_phys => kind_phot
  use module_prates_tuv, only: calc_tuv_init, calc_tuv_prates
  use params_mod, only: input_data_root

  implicit none

  integer, protected :: tuv_n_wavelen = -1

  integer, protected :: tuv_n_phot

  character(len=512) :: xsqy_filepath = 'NOTSET'
  character(len=512) :: etf_filepath = 'NOTSET'
  logical :: is_full_tuv = .true.
  
contains

  subroutine tuv_photolysis_readnl(nml_file, errmsg, errflg) ! this will be a CPF interface someday

    use module_prates_tuv, only: read_etf
    use wavelength_grid, only: nwave
    use wavelength_grid, only: wavelength_grid_init
 
    character(len=*), intent(in)  :: nml_file
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=512) :: wavelen_grid_filepath
    character(len=64) :: xsqy_file = 'NONE'
    character(len=64) :: etf_file = 'NONE'
    character(len=64) :: wavelength_grid_file = 'NONE'

    namelist /tuv_opts/ input_data_root, wavelength_grid_file, etf_file, is_full_tuv, xsqy_file

    open(unit=10,file=nml_file)
    read(unit=10,nml=tuv_opts)
    close(10)

    if (etf_file=='NONE') then
       etf_file = trim(wavelength_grid_file)
    end if

    wavelen_grid_filepath = trim(input_data_root)//'/'//trim(wavelength_grid_file)
    etf_filepath = trim(input_data_root)//'/'//trim(etf_file)
    xsqy_filepath = trim(input_data_root)//'/'//trim(xsqy_file)
    
    call wavelength_grid_init(wavelen_grid_filepath, errmsg, errflg)
    if(errflg/=0) then
       return
    end if

    call read_etf(etf_filepath, errmsg, errflg)
    if(errflg/=0) then
       return
    end if

    tuv_n_wavelen = nwave
    
  end subroutine tuv_photolysis_readnl

!> \section arg_table_tuv_photolysis_init Argument Table
!! \htmlinclude tuv_photolysis_init.html
!!
subroutine tuv_photolysis_init( realkind, nlev, jnames, tuv_n_wavelen, errmsg, errflg )
    use wavelength_grid, only: nwave

    integer,          intent(in)  :: realkind
    integer,          intent(in)  :: nlev
    character(len=*), intent(in)  :: jnames(:)

    !! NOTE THIS VARIABLE WILL GO AWAY - FOR NOW IS REQUIRED WORKAROUND FOR CPF
    integer,          intent(out) :: tuv_n_wavelen

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    if ( realkind/=kind_phys ) then
       errmsg = 'tuv_photolysis_init: realkind does not match kind_phot'
       errflg = 1
       return
    end if

    tuv_n_wavelen = nwave
    tuv_n_phot = size(jnames)
 
    call  calc_tuv_init( is_full_tuv, nlev, jnames, xsqy_filepath, errmsg, errflg )

  end subroutine tuv_photolysis_init

!> \section arg_table_tuv_photolysis_run Argument Table
!! \htmlinclude tuv_photolysis_run.html
!!
subroutine tuv_photolysis_run( nlev, temp, press_mid, radfld, srb_o2_xs, tuv_prates, errmsg, errflg )

    integer,          intent(in)  :: nlev
    real(kind_phys),  intent(in)  :: temp(:)
    real(kind_phys),  intent(in)  :: press_mid(:)
    real(kind_phys),  intent(in)  :: radfld(:,:) ! (nwave,nlev)
    real(kind_phys),  intent(in)  :: srb_o2_xs(:,:) !(nwave,kts:kte)
    real(kind_phys),  intent(out) :: tuv_prates(:,:) ! /sec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k, kk, j
    real(kind_phys) :: airdens(nlev) ! # molecules / cm3 in each layer
    real(kind_phys) :: tlev(nlev) ! # K -- bottom up

    real(kind_phys), parameter :: kboltz= 1.38064852e-16_kind_phys ! boltzmann constant (erg/K)

    ! inputs need to be bottom vertical coord
    do k=1,nlev
       kk=nlev-k+1
       airdens(kk) = 10._kind_phys*press_mid(k)/(kboltz*temp(k))
    end do
    tlev(nlev:1:-1) = temp(1:nlev)

    call calc_tuv_prates(1 ,nlev,nlev, tlev, airdens, radfld, srb_o2_xs, tuv_prates, errmsg, errflg)

    ! return top down rates
    do j=1,tuv_n_phot
       tuv_prates(:nlev,j) = tuv_prates(nlev:1:-1,j)
    end do
    
  end subroutine tuv_photolysis_run

end module tuv_photolysis
