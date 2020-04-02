module tuv_photolysis

  use environ_conditions_mod, only: environ_conditions_create, environ_conditions
  use phot_kind_mod, only: kind_phys => kind_phot
  use module_prates_tuv, only: calc_tuv_init, calc_tuv_prates
  use params_mod, only: input_data_root

  implicit none

  integer, protected :: tuv_n_wavelen = -1

  integer, protected :: tuv_n_phot

  character(len=512) :: xsqy_filepath = 'NOTSET'
  character(len=512) :: etf_filepath = 'NOTSET'
  logical :: is_full_tuv = .true.
  logical :: read_rate_constants_from_file = .false.

  !> Photolysis rates (s-1) dimesions: (number of boxes)
  type(environ_conditions),allocatable :: photo_rate_constants(:)

contains

  subroutine tuv_photolysis_readnl(nml_file, nbox, jnames, env_lat, env_lon, env_lev, errmsg, errflg) ! this will be a CPF interface someday

    use module_prates_tuv, only: read_etf
    use wavelength_grid, only: nwave
    use wavelength_grid, only: wavelength_grid_init
 
    character(len=*), intent(in)  :: nml_file
    integer         , intent(in)  :: nbox
    character(len=*), intent(in)  :: jnames(:)
    real            , intent(in)  :: env_lat(:)
    real            , intent(in)  :: env_lon(:)
    real            , intent(in)  :: env_lev(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! If a rate constants file is included in the tuv options, rate constants will be
    ! read from the file, otherwise they will be calculated
    character(len=512) :: rate_constants_file = '' ! Calculate rate constants by default
    character(len=512) :: wavelen_grid_filepath
    character(len=64) :: xsqy_file = 'NONE'
    character(len=64) :: etf_file = 'NONE'
    character(len=64) :: wavelength_grid_file = 'NONE'

    integer :: ibox

    namelist /tuv_opts/ rate_constants_file, input_data_root, wavelength_grid_file, etf_file, is_full_tuv, xsqy_file

    ! TUV setup
    open(unit=10,file=nml_file)
    read(unit=10,nml=tuv_opts)
    close(10)

    ! Load photolysis rate constants from a file, if indicated
    if (len(trim(rate_constants_file)).gt.0) then
      read_rate_constants_from_file = .true.
      allocate(photo_rate_constants(nbox))
      do ibox = 1, nbox
        photo_rate_constants(ibox) = environ_conditions_create( rate_constants_file,   &
                                                                lat = env_lat( ibox ), &
                                                                lon = env_lon( ibox ), &
                                                                lev = env_lev( ibox ) )
      end do
      write(*,*) "Reading photolysis rates from file '"//trim(rate_constants_file)//"'"
    end if

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
subroutine tuv_photolysis_init( realkind, nlyr, jnames, tuv_n_wavelen, errmsg, errflg )
    use wavelength_grid, only: nwave

    integer,          intent(in)  :: realkind
    integer,          intent(in)  :: nlyr
    character(len=*), intent(in)  :: jnames(:)
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
 
    call  calc_tuv_init( is_full_tuv, nlyr, jnames, xsqy_filepath, errmsg, errflg )

  end subroutine tuv_photolysis_init

!> \section arg_table_tuv_photolysis_run Argument Table
!! \htmlinclude tuv_photolysis_run.html
!!
subroutine tuv_photolysis_run( TimeStart, nlyr, nbox, jnames, temp, press_mid, radfld, srb_o2_xs, tuv_prates, errmsg, errflg )

    real(kind_phys),  intent(in)  :: TimeStart
    integer,          intent(in)  :: nlyr
    integer,          intent(in)  :: nbox
    character(len=*), intent(in)  :: jnames(:)
    real(kind_phys),  intent(in)  :: temp(:)
    real(kind_phys),  intent(in)  :: press_mid(:)
    real(kind_phys),  intent(in)  :: radfld(:,:) ! (nwave,nlyr)
    real(kind_phys),  intent(in)  :: srb_o2_xs(:,:) !(nwave,kts:kte)
    real(kind_phys),  intent(out) :: tuv_prates(:,:) ! /sec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: k, kk, j, i_photo_rxn
    real(kind_phys) :: airdens(nlyr) ! # molecules / cm3 in each layer
    real(kind_phys) :: tlev(nlyr) ! # K -- bottom up

    real(kind_phys), parameter :: kboltz= 1.38064852e-16_kind_phys ! boltzmann constant (erg/K)

    ! use the photo rate constants from the input file, if indicated
    ! \todo the tuv_prates seems to be set up for only a single box
    if (read_rate_constants_from_file) then
      call photo_rate_constants(1)%update(TimeStart)
      do i_photo_rxn = 1, size(jnames)
        tuv_prates(:nlyr, i_photo_rxn) = &
          photo_rate_constants(1)%getvar("photo_rate_constant_"//trim(jnames(i_photo_rxn)))
      end do
      return
    end if

    ! inputs need to be bottom vertical coord
    do k=1,nlyr
       kk=nlyr-k+1
       airdens(kk) = 10._kind_phys*press_mid(k)/(kboltz*temp(k))
    end do
    tlev(nlyr:1:-1) = temp(1:nlyr)

    call calc_tuv_prates(1 ,nlyr,nlyr, tlev, airdens, radfld, srb_o2_xs, tuv_prates, errmsg, errflg)

    ! return top down rates
    do j=1,tuv_n_phot
       tuv_prates(:nlyr,j) = tuv_prates(nlyr:1:-1,j)
    end do
    
  end subroutine tuv_photolysis_run

end module tuv_photolysis
