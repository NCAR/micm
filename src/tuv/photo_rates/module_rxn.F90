!=============================================================================*
! This file contains the following subroutines, related to reading/loading
! the product (cross section) x (quantum yield) for photo-reactions:
!     r01 through r47
!     r101 through r148, skipped r116,r117, added pxCH2O
!=============================================================================*
      module module_rxn

      use phot_kind_mod, only: rk => kind_phot
      use params_mod, only: largest, deltax, input_data_root, kin, qnan
      use  numer_mod, only: addpnt, inter2, inter4
      
      implicit none

      private :: fo3qy2, qyacet

      logical, private :: initialize = .true.

      ! internal array dimensions
      integer, parameter :: max_files = 5
      ! altitude, wavelength, time (or solar zenith angle) grids
      ! altitude
      integer, parameter :: kz=125
      ! wavelength
      integer, parameter :: kw=1000
      !  wavelength and altitude dependent
      integer, parameter :: kj=150

      ! small numbers (positive and negative)
      real(rk), parameter :: pzero = +10._rk/largest

      integer :: npht, npht_tab
      
      type file_specs
        integer            :: nfiles
        integer            :: nskip(max_files)
        integer            :: nread(max_files)
        real(rk)           :: xfac(max_files)
        character(len=388) :: filename(max_files)
      end type file_specs

      type xs_qy_tab
        integer :: tpflag
        integer :: channel
        integer :: jndx
        real(rk)    :: qyld
        character(len=50) :: equation
        character(len=50) :: rxn_name
        type(xs_qy_tab), pointer :: next
        type(xs_qy_tab), pointer :: last
        type(file_specs)  :: filespec
      end type xs_qy_tab

      type(xs_qy_tab), allocatable, target :: xsqy_tab(:)
      type(xs_qy_tab), pointer             :: xsqy_tab_head
      type(xs_qy_tab), pointer             :: xsqy_tab_tail

!=====================================================================
!  the following is fortran2003 compliant code
!=====================================================================
      type xsqy_subs
        procedure(xsqy), nopass, pointer :: xsqy_sub
      end type xsqy_subs

      abstract interface
        SUBROUTINE xsqy(nw,wl,wc,nz,tlev,airden,j,errmsg,errflg, sq)

          use phot_kind_mod, only: rk => kind_phot
          
          INTEGER, intent(in) :: nw
          INTEGER, intent(in) :: nz
          REAL(rk), intent(in)    :: wl(:), wc(:)
          REAL(rk), intent(in)    :: tlev(:)
          REAL(rk), intent(in)    :: airden(:)
          character(len=*), intent(out) :: errmsg
          integer,          intent(out) :: errflg
          real(rk), optional, intent(out) :: sq(:,:)
          
          INTEGER, intent(inout) :: j
        end SUBROUTINE xsqy
      end interface

      type(xsqy_subs), allocatable :: the_subs(:)
      real(rk) :: xnan

      real(rk), parameter :: wlla(2)  = (/ 121.4_rk, 121.9_rk/) ! Lyamn Alpha limits
      integer :: la_ndx = -1
      
      CONTAINS

      SUBROUTINE no_z_dep(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!  generic routine
!-----------------------------------------------------------------------------*

!      use module_params

! input

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk),    intent(in) :: wl(:), wc(:)
      REAL(rk),    intent(in) :: airden(:)
      REAL(rk),    intent(in) :: tlev(:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)
      
      integer, PARAMETER :: kdata=600
      integer, PARAMETER :: jdata=200

! local
      REAL(rk) :: x1(kdata)
      REAL(rk) :: y1(kdata)
      REAL(rk) :: yg(kw)
      real(rk), save :: xsq(kw,jdata)

      errmsg = ' '
      errflg = 0

      if (present(sq)) then
         sq = xnan
      end if

      if( initialize ) then
        yg = xnan
        CALL readit
        if( xsqy_tab(j)%qyld == 1._rk ) then
!*** quantum yield assumed to be unity
          xsq(1:nw-1,j) = yg(1:nw-1)
        else
          xsq(1:nw-1,j) = xsqy_tab(j)%qyld * yg(1:nw-1)
        endif
      else
         sq(1:nw-1,1) = xsq(1:nw-1,j)
      endif

      CONTAINS 

      SUBROUTINE readit

      integer :: n, fileno
      character(len=132) :: filename

      do fileno = 1,xsqy_tab(j)%filespec%nfiles
        filename = trim( xsqy_tab(j)%filespec%filename(fileno) )
        n = xsqy_tab(j)%filespec%nread(fileno)
        if( xsqy_tab(j)%filespec%nskip(fileno) >= 0 ) then
          CALL base_read( filespec=trim(filename), errmsg=errmsg, errflg=errflg, &
                          skip_cnt=xsqy_tab(j)%filespec%nskip(fileno), &
                          rd_cnt  =n,x=x1,y=y1 )
        else
          CALL base_read( filespec=trim(filename), errmsg=errmsg, errflg=errflg, rd_cnt=n,x=x1,y=y1 )
        endif
        y1(1:n) = y1(1:n) * xsqy_tab(j)%filespec%xfac(fileno)

        CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      enddo

      END SUBROUTINE readit

      END SUBROUTINE no_z_dep

      LOGICAL FUNCTION get_initialization()

      get_initialization = initialize

      END FUNCTION get_initialization

      SUBROUTINE set_initialization( status )

      LOGICAL, intent(in) :: status

      initialize = status

      END SUBROUTINE set_initialization

      SUBROUTINE rxn_init( nw, wl, errmsg, errflg )
!---------------------------------------------
!  initialize wrf-tuv
!---------------------------------------------

      integer, intent(in) :: nw
      real(rk), intent(in)    :: wl(nw)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: astat, m, iw
      real(rk) :: wc(nw-1)
      
      errmsg = ' '
      errflg = 0
      xnan = qnan()

      ! find lyman alpha band
      wc(1:nw-1) = 0.5_rk*(wl(1:nw-1)+wl(2:nw))
      do iw = 1, nw
         if (wc(iw)>wlla(1) .and. wc(iw)<wlla(2) .and. wc(iw+1)>wlla(2)) then
            la_ndx = iw
            exit
         end if
      end do
      if (la_ndx < 1) then
         errflg = 1
         errmsg = 'rxn_init: Not able to find Lyamn Alpha wavelength band'
         return
      end if

      call set_initialization( status=.true. )

      if( .not. allocated( xsqy_tab ) ) then
         allocate( xsqy_tab(kj),stat=astat )
         if( astat /= 0 ) then
            write(errmsg,'(''rxn_init: failed to allocate xsqy_tab; error = '',i4)') astat
            errflg = astat
            return
         endif
      endif
      if( .not. allocated( the_subs ) ) then
         allocate( the_subs(kj),stat=astat )
         if( astat /= 0 ) then
            write(errmsg,'(''rxn_init: failed to allocate xsqy_tab subs; error = '',i4)') astat
            errflg = astat
            return
         endif
      endif

      nullify( xsqy_tab_head )
      nullify( xsqy_tab_tail )
      
      xsqy_tab(1:kj)%tpflag  = 0
      xsqy_tab(1:kj)%channel = 1
      xsqy_tab(1:kj)%equation =  ' '
      xsqy_tab(1:kj)%qyld    =  1._rk
      xsqy_tab(1:kj)%filespec%nfiles =  1
      do m = 1,max_files
        xsqy_tab(1:kj)%filespec%nskip(m) =  0
        xsqy_tab(1:kj)%filespec%nread(m) =  0
        xsqy_tab(1:kj)%filespec%xfac(m)  =  1.e-20_rk
        xsqy_tab(1:kj)%filespec%filename(m) = ' '
      end do
      do m = 1,kj
        nullify( xsqy_tab(m)%next )
        nullify( xsqy_tab(m)%last )
        the_subs(m)%xsqy_sub => null()
      end do

      npht_tab = 2
      call setup_sub_calls( the_subs,npht_tab )

      call diagnostics

      END SUBROUTINE rxn_init

      subroutine setup_sub_calls( subr, m )

      integer, intent(inout) :: m
      type(xsqy_subs), intent(inout) :: subr(:)

      xsqy_tab(m)%equation   = 'O2 + hv -> O(1D) + O(3P)'
      xsqy_tab(m)%rxn_name   = 'jo2_a'
      xsqy_tab(m+1)%equation = 'O2 + hv -> O(3P) + O(3P)'
      xsqy_tab(m+1)%rxn_name = 'jo2_b'
      xsqy_tab(m:m+1)%jndx = (/ m,m+1 /)
      xsqy_tab(m:m+1)%channel = (/1,2/)
      xsqy_tab(m:m+1)%tpflag = 0 ! quantum yield is not temperature nor pressure dependent 
      subr(m  )%xsqy_sub => qy_o2
      subr(m+1)%xsqy_sub => qy_o2
      m = m + 2

      xsqy_tab(m)%equation   = 'O3 -> O2 + O(1D)'
      xsqy_tab(m)%rxn_name   = 'j_o3_a'
      xsqy_tab(m+1)%equation = 'O3 -> O2 + O(3P)'
      xsqy_tab(m+1)%rxn_name = 'j_o3_b'
      xsqy_tab(m:m+1)%jndx = (/ m,m+1 /)
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 1
      subr(m  )%xsqy_sub => r01
      subr(m+1)%xsqy_sub => r01
      m = m + 2

      xsqy_tab(m)%equation   = 'O3 -> O2 + O(1D)'
      xsqy_tab(m)%rxn_name   = 'jo3_a'
      xsqy_tab(m+1)%equation = 'O3 -> O2 + O(3P)'
      xsqy_tab(m+1)%rxn_name = 'jo3_b'
      xsqy_tab(m:m+1)%jndx = (/ m,m+1 /)
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 1
      xsqy_tab(m:m+1)%filespec%nfiles = 3
      xsqy_tab(m:m+1)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_O3_Ackerman_1971.txt'
      xsqy_tab(m:m+1)%filespec%nskip(1) = 9
      xsqy_tab(m:m+1)%filespec%nread(1) = 56
      xsqy_tab(m:m+1)%filespec%filename(2) = trim(input_data_root)//'/XSQY/XS_O3_218_JPL06.txt'
      xsqy_tab(m:m+1)%filespec%nskip(2) = 28
      xsqy_tab(m:m+1)%filespec%nread(2) = 65
      xsqy_tab(m:m+1)%filespec%filename(3) = trim(input_data_root)//'/XSQY/XS_O3_298_JPL06.txt'
      xsqy_tab(m:m+1)%filespec%nskip(3) = 39
      xsqy_tab(m:m+1)%filespec%nread(3) = 168
      subr(m  )%xsqy_sub => XSQY_O3
      subr(m+1)%xsqy_sub => XSQY_O3
      m = m + 2

      xsqy_tab(m)%equation = 'NO2 -> NO + O(3P)'
      xsqy_tab(m)%rxn_name = 'jno2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/YLD/NO2_jpl11.yld'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 25
      subr(m)%xsqy_sub   => r02
      m = m + 1

      xsqy_tab(m)%equation   = 'NO3 -> NO + O2'
      xsqy_tab(m)%rxn_name   = 'j_no3_b'
      xsqy_tab(m+1)%equation = 'NO3 -> NO2 + O(3P)'
      xsqy_tab(m+1)%rxn_name = 'j_no3_a'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 1
      xsqy_tab(m:m+1)%filespec%nfiles = 2
      xsqy_tab(m:m+1)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/NO3_jpl11.abs'
      xsqy_tab(m:m+1)%filespec%nskip(1) = 6
      xsqy_tab(m:m+1)%filespec%nread(1) = 289
      xsqy_tab(m:m+1)%filespec%filename(2) = trim(input_data_root)//'/DATAJ1/YLD/NO3_jpl2011.qy'
      xsqy_tab(m:m+1)%filespec%nskip(2) = 5
      xsqy_tab(m:m+1)%filespec%nread(2) = 56
      xsqy_tab(m:m+1)%filespec%xfac(2)  = 1.e-3_rk
      subr(m)%xsqy_sub   => r03
      subr(m+1)%xsqy_sub => r03
      m = m + 2

      
      xsqy_tab(m)%equation   = 'NO3 -> NO + O2'
      xsqy_tab(m)%rxn_name   = 'jno3_b'
      xsqy_tab(m+1)%equation = 'NO3 -> NO2 + O(3P)'
      xsqy_tab(m+1)%rxn_name = 'jno3_a'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 1
      xsqy_tab(m:m+1)%filespec%nfiles = 2
      xsqy_tab(m:m+1)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_NO3_JPL06.txt'
      xsqy_tab(m:m+1)%filespec%nskip(1) = 26
      xsqy_tab(m:m+1)%filespec%nread(1) = 289
      xsqy_tab(m:m+1)%filespec%filename(2) = trim(input_data_root)//'/XSQY/QY_NO3_JPL06.txt'
      xsqy_tab(m:m+1)%filespec%nskip(2) = 23
      xsqy_tab(m:m+1)%filespec%nread(2) = 56
      xsqy_tab(m:m+1)%filespec%xfac(2)  = 1.e-3_rk
      subr(m)%xsqy_sub   => XSQY_NO3
      subr(m+1)%xsqy_sub => XSQY_NO3
      m = m + 2

      xsqy_tab(m)%equation   = 'N2O5 -> NO3 + NO + O(3P)'
      xsqy_tab(m)%rxn_name   = 'jn2o5_b'
      xsqy_tab(m+1)%equation = 'N2O5 -> NO3 + NO2'
      xsqy_tab(m+1)%rxn_name = 'jn2o5_a'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 1
      xsqy_tab(m)%filespec%nfiles = 1
      xsqy_tab(m:m+1)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_N2O5_JPL06.txt'
      xsqy_tab(m:m+1)%filespec%nskip(1) = 41
      xsqy_tab(m:m+1)%filespec%nread(1) = 103
      subr(m)%xsqy_sub   => XSQY_N2O5
      subr(m+1)%xsqy_sub => XSQY_N2O5
      m = m + 2

      xsqy_tab(m  )%equation   = 'H2O + hv -> H + OH'
      xsqy_tab(m  )%rxn_name   = 'jh2o_a'
      xsqy_tab(m+1)%equation   = 'H2O + hv -> H2 + O(1D)'
      xsqy_tab(m+1)%rxn_name   = 'jh2o_b'
      xsqy_tab(m+2)%equation   = 'H2O + hv -> 2H + O(3P)'
      xsqy_tab(m+2)%rxn_name   = 'jh2o_c'
      xsqy_tab(m  )%jndx = m
      xsqy_tab(m+1)%jndx = m+1
      xsqy_tab(m+2)%jndx = m+2
      xsqy_tab(m  )%channel = 1
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m+2)%channel = 3
      xsqy_tab(m:m+2)%tpflag = 0
      xsqy_tab(m:m+2)%filespec%nfiles = 3
      xsqy_tab(m:m+2)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_H2O_JPL06.txt'
      xsqy_tab(m:m+2)%filespec%nskip(1) = 28
      xsqy_tab(m:m+2)%filespec%nread(1) = 8
      xsqy_tab(m:m+2)%filespec%filename(2) = trim(input_data_root)//'/XSQY/XS_H2O_cantrell_1996.txt'
      xsqy_tab(m:m+2)%filespec%nskip(2) = 12
      xsqy_tab(m:m+2)%filespec%nread(2) = 11
      xsqy_tab(m:m+2)%filespec%filename(3) = trim(input_data_root)//'/XSQY/XS_H2O_yoshino_1996.txt'
      xsqy_tab(m:m+2)%filespec%nskip(3) = 23
      xsqy_tab(m:m+2)%filespec%nread(3) = 6783
      xsqy_tab(m  )%filespec%xfac(1:3) = (/ 1.e-20_rk, 1.e-20_rk, 1._rk /)
      xsqy_tab(m+1)%filespec%xfac(1:3) = (/ 1.e-20_rk, 1.e-20_rk, 1._rk /)
      xsqy_tab(m+2)%filespec%xfac(1:3) = (/ 1.e-20_rk, 1.e-20_rk, 1._rk /)
      subr(m  )%xsqy_sub => XSQY_H2O
      subr(m+1)%xsqy_sub => XSQY_H2O
      subr(m+2)%xsqy_sub => XSQY_H2O
      m = m + 3

      xsqy_tab(m  )%channel = 1
      xsqy_tab(m  )%equation   = 'HO2NO2 -> OH + NO3'
      xsqy_tab(m  )%rxn_name   = 'jho2no2_a'
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m+1)%equation   = 'HO2NO2 -> HO2 + NO2'
      xsqy_tab(m+1)%rxn_name   = 'jho2no2_b'
      xsqy_tab(m  )%jndx = m
      xsqy_tab(m+1)%jndx = m+1
      xsqy_tab(m:m+2)%tpflag = 1
      xsqy_tab(m:m+2)%filespec%nfiles = 2
      xsqy_tab(m:m+2)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_HO2NO2_JPL06.txt'
      xsqy_tab(m:m+2)%filespec%nskip(1) = 30
      xsqy_tab(m:m+2)%filespec%nread(1) = 54     
      xsqy_tab(m:m+2)%filespec%filename(2) = trim(input_data_root)//'/XSQY/XS_HO2NO2_JPL06.txt'
      xsqy_tab(m:m+2)%filespec%nskip(2) = 86
      xsqy_tab(m:m+2)%filespec%nread(2) = 36
      subr(m  )%xsqy_sub => XSQY_HO2NO2
      subr(m+1)%xsqy_sub => XSQY_HO2NO2
      m = m + 2

      xsqy_tab(m)%jndx = m
      xsqy_tab(m)%equation   = 'CH4 + hv -> CH3O2 + H'
      xsqy_tab(m)%rxn_name   = 'jch4_a'
      xsqy_tab(m)%channel = 1
      xsqy_tab(m+1)%jndx = m+1
      xsqy_tab(m+1)%equation = 'CH4 + hv -> H2 + Products'
      xsqy_tab(m+1)%rxn_name = 'jch4_b'
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 0
      xsqy_tab(m:m+1)%filespec%nfiles = 1
      xsqy_tab(m:m+1)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_CH4.txt'
      xsqy_tab(m:m+1)%filespec%nskip(1) = 22
      xsqy_tab(m:m+1)%filespec%nread(1) = 39
      xsqy_tab(m:m+1)%filespec%xfac(1)  = 1._rk
      xsqy_tab(m  )%qyld  = 0.45_rk
      xsqy_tab(m+1)%qyld  = 0.55_rk
      subr(m)%xsqy_sub   => no_z_dep
      subr(m+1)%xsqy_sub => no_z_dep
      m = m + 2

      xsqy_tab(m)%equation = 'CO2 + hv -> CO + O'
      xsqy_tab(m)%rxn_name = 'jco2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_CO2.txt'
      xsqy_tab(m)%filespec%nskip(1) = 20
      xsqy_tab(m)%filespec%nread(1) = 55
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'COF2 + hv  -> 2F'
      xsqy_tab(m)%rxn_name = 'jcof2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_COF2_JPL10.txt'
      xsqy_tab(m)%filespec%nskip(1) = 31
      xsqy_tab(m)%filespec%nread(1) = 21
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'COFCl + hv  -> Cl + F'
      xsqy_tab(m)%rxn_name = 'jcofcl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_COFCL_JPL10.txt'
      xsqy_tab(m)%filespec%nskip(1) = 32
      xsqy_tab(m)%filespec%nread(1) = 32
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'HBr + hv -> H + Br'
      xsqy_tab(m)%rxn_name = 'jhbr'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_HBR_JPL06.txt'
      xsqy_tab(m)%filespec%nskip(1) = 44
      xsqy_tab(m)%filespec%nread(1) = 40
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'HF + hv  -> H + F'
      xsqy_tab(m)%rxn_name = 'jhf'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_HF.txt'
      xsqy_tab(m)%filespec%nskip(1) = 14
      xsqy_tab(m)%filespec%nread(1) = 39
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'SF6 + hv -> product'
      xsqy_tab(m)%rxn_name = 'jsf6'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_SF6.txt'
      xsqy_tab(m)%filespec%nskip(1) = 14
      xsqy_tab(m)%filespec%nread(1) = 14
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'GLYALD + hv -> 2HO2 + CO + CH2O' !  GLYALD (HOCH2CHO) + hv -> 2*HO2 + CO + CH2O 
      xsqy_tab(m)%rxn_name = 'jglyald'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m)%qyld = 0.5_rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_GLYALD.txt'
      xsqy_tab(m)%filespec%nskip(1) = 15
      xsqy_tab(m)%filespec%nread(1) = 131
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'HYAC + hv -> CH3CO3 + HO2 + CH2O'
      xsqy_tab(m)%rxn_name = 'jhyac'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m)%qyld = 0.65_rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_HYAC.txt'
      xsqy_tab(m)%filespec%nskip(1) = 8
      xsqy_tab(m)%filespec%nread(1) = 101
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'MACR + hv -> 1.34HO2 + 0.66MCO3 + 1.34CH2O+ CH3CO3'
      xsqy_tab(m)%rxn_name = 'jmacr_a'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m)%qyld = 0.005_rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_MACR_JPL06.txt'
      xsqy_tab(m)%filespec%nskip(1) = 20
      xsqy_tab(m)%filespec%nread(1) = 146
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'MACR + hv -> 0.66OH + 1.34CO'
      xsqy_tab(m)%rxn_name = 'jmacr_b'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m)%qyld = 0.005_rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_MACR_JPL06.txt'
      xsqy_tab(m)%filespec%nskip(1) = 20
      xsqy_tab(m)%filespec%nread(1) = 146
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'H2SO4 + hv -> HSO3 + OH'
      xsqy_tab(m)%rxn_name = 'jh2so4'
      xsqy_tab(m)%jndx   = m
      xsqy_tab(m)%tpflag = 0
      xsqy_tab(m)%qyld   = 1._rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_H2SO4_mills.txt'
      xsqy_tab(m)%filespec%nskip(1) = 13
      xsqy_tab(m)%filespec%nread(1) = 125
      subr(m)%xsqy_sub => XSQY_MMILLS
      m = m + 1

      xsqy_tab(m)%equation = 'OCS + hv -> CO + S'
      xsqy_tab(m)%rxn_name = 'jocs'
      xsqy_tab(m)%jndx   = m
      xsqy_tab(m)%tpflag = 0
      xsqy_tab(m)%qyld   = 1._rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_OCS_mills.txt'
      xsqy_tab(m)%filespec%nskip(1) = 13
      xsqy_tab(m)%filespec%nread(1) = 125
      subr(m)%xsqy_sub => XSQY_MMILLS
      m = m + 1

      xsqy_tab(m)%equation = 'SO + hv -> S + O'
      xsqy_tab(m)%rxn_name = 'jso'
      xsqy_tab(m)%jndx   = m
      xsqy_tab(m)%tpflag = 0
      xsqy_tab(m)%qyld   = 1._rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_SO_mills.txt'
      xsqy_tab(m)%filespec%nskip(1) = 13
      xsqy_tab(m)%filespec%nread(1) = 125
      subr(m)%xsqy_sub => XSQY_MMILLS
      m = m + 1

      xsqy_tab(m)%equation = 'SO2 + hv -> SO + O'
      xsqy_tab(m)%rxn_name = 'jso2'
      xsqy_tab(m)%jndx   = m
      xsqy_tab(m)%tpflag = 0
      xsqy_tab(m)%qyld   = 1._rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_SO2_mills.txt'
      xsqy_tab(m)%filespec%nskip(1) = 13
      xsqy_tab(m)%filespec%nread(1) = 125
      subr(m)%xsqy_sub => XSQY_MMILLS
      m = m + 1

      xsqy_tab(m)%equation = 'SO3 + hv -> SO2 + O'
      xsqy_tab(m)%rxn_name = 'jso3'
      xsqy_tab(m)%jndx   = m
      xsqy_tab(m)%tpflag = 0
      xsqy_tab(m)%qyld   = 1._rk
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_SO3_mills.txt'
      xsqy_tab(m)%filespec%nskip(1) = 13
      xsqy_tab(m)%filespec%nread(1) = 125
      subr(m)%xsqy_sub => XSQY_MMILLS
      m = m + 1

      xsqy_tab(m)%equation = 'CH2Br2 + hv -> 2Br'
      xsqy_tab(m)%rxn_name = 'jch2br2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_CH2BR2_JPL06.txt'
      xsqy_tab(m)%filespec%nskip(1) = 53
      xsqy_tab(m)%filespec%nread(1) = 34
      subr(m)%xsqy_sub   => XSQY_CH2BR2
      m = m + 1

      xsqy_tab(m)%equation = 'NO = hv => NOp + e'
      xsqy_tab(m)%rxn_name = 'jno_i'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 0
      subr(m  )%xsqy_sub => XSQY_NOp
      m = m + 1

      xsqy_tab(m)%equation = 'HNO2 -> OH + NO'
      xsqy_tab(m)%rxn_name = 'jhno2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HONO_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 192
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'HNO3 -> OH + NO2'
      xsqy_tab(m)%rxn_name = 'jhno3'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HNO3_burk.abs'
      xsqy_tab(m)%filespec%nskip(1) = 6
      xsqy_tab(m)%filespec%nread(1) = 83
      subr(m)%xsqy_sub   => r06
      m = m + 1

      xsqy_tab(m)%equation = 'HNO4 -> HO2 + NO2'
      xsqy_tab(m)%rxn_name = 'jhno4'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HNO4_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 54
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'H2O2 -> 2 OH'
      xsqy_tab(m)%rxn_name = 'jh2o2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%nfiles      = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/H2O2_jpl94.abs'
      xsqy_tab(m)%filespec%filename(2) = trim(input_data_root)//'/DATAJ1/ABS/H2O2_Kahan.abs'
      xsqy_tab(m)%filespec%nskip(1:2) = (/ -1,0 /)
      xsqy_tab(m)%filespec%nread(2)   = 494
      subr(m)%xsqy_sub   => r08
      m = m + 1

      xsqy_tab(m)%equation = 'CHBr3 -> Products'
      xsqy_tab(m)%rxn_name = 'jchbr3'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CHBr3.jpl97'
      xsqy_tab(m)%filespec%nskip(1) = 6
      xsqy_tab(m)%filespec%nread(1) = 87
      subr(m)%xsqy_sub   => r09
      m = m + 1

      xsqy_tab(m)%equation = 'CH3CHO + hv -> CH3O2 + CO + HO2'
      xsqy_tab(m)%rxn_name = 'jch3cho'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 2
      xsqy_tab(m)%filespec%nfiles = 3
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_CH3CHO_JPL06.txt'
      xsqy_tab(m)%filespec%filename(2) = trim(input_data_root)//'/XSQY/QY_CH3CHO_iup.yld.txt'
      xsqy_tab(m)%filespec%filename(3) = trim(input_data_root)//'/XSQY/QY_CH3CHO_press.yld.txt'
      xsqy_tab(m)%filespec%nskip(1) = 31
      xsqy_tab(m)%filespec%nread(1) = 101
      xsqy_tab(m)%filespec%nskip(2) = 4
      xsqy_tab(m)%filespec%nread(2) = 12
      xsqy_tab(m)%filespec%nskip(3) = 4
      xsqy_tab(m)%filespec%nread(3) = 5
      subr(m)%xsqy_sub   => XSQY_CH3CHO
      m = m + 1

      xsqy_tab(m)%equation   = 'CH3CHO -> CH3 + HCO'
      xsqy_tab(m+1)%equation = 'CH3CHO -> CH4 + CO'
      xsqy_tab(m+2)%equation = 'CH3CHO -> CH3CO + H'
      xsqy_tab(m)%rxn_name = 'j_ch3cho_a'
      xsqy_tab(m+1)%rxn_name = 'j_ch3cho_b'
      xsqy_tab(m+2)%rxn_name = 'j_ch3cho_c'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1:m+2)%channel = (/ 2,3 /)
      xsqy_tab(m:m+2)%tpflag = (/ 2,0,0 /)
      xsqy_tab(m)%filespec%nfiles = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CH3CHO/CH3CHO_jpl11.abs'
      xsqy_tab(m)%filespec%filename(2) = trim(input_data_root)//'/DATAJ1/CH3CHO/CH3CHO_uip.yld'
      xsqy_tab(m)%filespec%nskip(1:2) = (/ 2,4 /)
      xsqy_tab(m)%filespec%nread(1:2) = (/ 101,12 /)
      subr(m)%xsqy_sub   => r11
      subr(m+1)%xsqy_sub => r11
      subr(m+2)%xsqy_sub => r11
      m = m + 3

      xsqy_tab(m)%equation = 'C2H5CHO -> C2H5 + HCO'
      xsqy_tab(m)%rxn_name = 'j_c2h5cho'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 2
      xsqy_tab(m)%filespec%nfiles = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/C2H5CHO/C2H5CHO_iup.abs'
      xsqy_tab(m)%filespec%filename(2) = trim(input_data_root)//'/DATAJ1/C2H5CHO/C2H5CHO_iup.yld'
      xsqy_tab(m)%filespec%nskip(1:2) = 4
      xsqy_tab(m)%filespec%nread(1:2) = (/ 106,5 /)
      subr(m)%xsqy_sub   => r12
      m = m + 1

      xsqy_tab(m)%equation   = 'CHOCHO -> HCO + HCO'
      xsqy_tab(m+1)%equation = 'CHOCHO -> H2 + 2CO'
      xsqy_tab(m+2)%equation = 'CHOCHO -> CH2O + CO'
      xsqy_tab(m)%rxn_name = 'j_gly_a'
      xsqy_tab(m+1)%rxn_name = 'j_gly_b'
      xsqy_tab(m+2)%rxn_name = 'j_gly_c'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1:m+2)%channel = (/ 2,3 /)
      xsqy_tab(m)%filespec%nfiles = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CHOCHO/glyoxal_jpl11.abs'
      xsqy_tab(m)%filespec%filename(2) = trim(input_data_root)//'/DATAJ1/CHOCHO/glyoxal_jpl11.qy'
      xsqy_tab(m)%filespec%nskip(1:2) = (/ 2,3 /)
      xsqy_tab(m)%filespec%nread(1:2) = (/ 277,40 /)
      subr(m)%xsqy_sub   => r13
      subr(m+1)%xsqy_sub => r13
      subr(m+2)%xsqy_sub => r13
      m = m + 3

      xsqy_tab(m)%equation = 'CH3COCHO -> CH3CO + HCO'
      xsqy_tab(m)%rxn_name = 'j_mgly'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CH3COCHO/CH3COCHO_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 294
      subr(m)%xsqy_sub   => r14
      m = m + 1

      xsqy_tab(m)%equation = 'MGLY + hv ->  CH3CO3  + CO + HO2'
      xsqy_tab(m)%rxn_name = 'jmgly'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_CH3COCHO_ncar.txt'
      xsqy_tab(m)%filespec%nskip(1) = 0
      xsqy_tab(m)%filespec%nread(1) = 271
      subr(m)%xsqy_sub   => XSQY_MGLY
      m = m + 1

      xsqy_tab(m)%equation = 'CH3COCH3 -> CH3CO + CH3'
      xsqy_tab(m)%rxn_name = 'j_ch3coch3'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 3
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CH3COCH3/CH3COCH3_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 5
      xsqy_tab(m)%filespec%nread(1) = 135
      subr(m)%xsqy_sub   => r15
      m = m + 1

      xsqy_tab(m)%equation = 'CH3COCH3 + hv -> CH3CO3 + CH3O2'
      xsqy_tab(m)%rxn_name = 'jacet'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%nfiles = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_ACETONE_JPL06.txt'
      xsqy_tab(m)%filespec%filename(2) = trim(input_data_root)//'/XSQY/XS_ACETONE_TD_JPL06.txt'
      xsqy_tab(m)%filespec%nskip(1:2) = (/ 35, 34 /)
      xsqy_tab(m)%filespec%nread(1:2) = (/ 135, 135 /)
      subr(m)%xsqy_sub   => XSQY_ACETONE
      m = m + 1

      xsqy_tab(m)%equation = 'CH3OOH -> CH3O + OH'
      xsqy_tab(m)%rxn_name = 'j_ch3ooh'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CH3OOH/CH3OOH_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 40
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3OOH + hv -> CH2O + H + OH'
      xsqy_tab(m)%rxn_name = 'jch3ooh'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_CH3OOH_JPL06.txt'
      xsqy_tab(m)%filespec%nskip(1) = 20
      xsqy_tab(m)%filespec%nread(1) = 32
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'Cl2O2 + hv -> Cl + ClOO'
      xsqy_tab(m)%rxn_name = 'jcl2o2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_CL2O2_JPL10_500nm.txt'
      xsqy_tab(m)%filespec%nskip(1) = 32
      xsqy_tab(m)%filespec%nread(1) = 521
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3ONO2 -> CH3O + NO2'
      xsqy_tab(m)%rxn_name = 'j_ch3ono2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/RONO2/CH3ONO2_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 65
      subr(m)%xsqy_sub   => r17
      m = m + 1

      xsqy_tab(m)%equation = 'PAN + hv -> Products'
      xsqy_tab(m)%rxn_name = 'jpan'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/XSQY/XS_PAN_JPL06.txt'
      xsqy_tab(m)%filespec%nskip(1) = 20
      xsqy_tab(m)%filespec%nread(1) = 78
      subr(m)%xsqy_sub   => XSQY_PAN
      m = m + 1

      xsqy_tab(m)%equation   = 'CH3CO(OONO2) -> CH3CO(OO) + NO2'
      xsqy_tab(m+1)%equation = 'CH3CO(OONO2) -> CH3CO(O) + NO3'
      xsqy_tab(m)%rxn_name = 'j_pan_a'
      xsqy_tab(m+1)%rxn_name = 'j_pan_b'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/RONO2/PAN_talukdar.abs'
      xsqy_tab(m)%filespec%nskip(1) = 14
      xsqy_tab(m)%filespec%nread(1) = 78
      subr(m)%xsqy_sub   => r18
      subr(m+1)%xsqy_sub => r18
      m = m + 2

      xsqy_tab(m)%equation = 'CCl2O -> Products'
      xsqy_tab(m)%rxn_name = 'j_ccl2o'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CCl2O_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CCl4 -> Products'
      xsqy_tab(m)%rxn_name = 'jccl4'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CCl4_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 5
      xsqy_tab(m)%filespec%nread(1) = 44
      subr(m)%xsqy_sub   => r20
      m = m + 1

      xsqy_tab(m)%equation = 'CClFO -> Products'
      xsqy_tab(m)%rxn_name = 'j_cclfo'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CClFO_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CF2O -> Products'
      xsqy_tab(m)%rxn_name = 'j_cf2o'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CF2O_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 5
      xsqy_tab(m)%filespec%nread(1) = 21
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CF2ClCFCl2 (CFC-113) -> Products'
      xsqy_tab(m)%rxn_name = 'jcfc113'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CFC-113_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => r23
      m = m + 1

      xsqy_tab(m)%equation = 'CF2ClCF2Cl (CFC-114) -> Products'
      xsqy_tab(m)%rxn_name = 'jcfc114'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CFC-114_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => r24
      m = m + 1

      xsqy_tab(m)%equation = 'CF3CF2Cl (CFC-115) -> Products'
      xsqy_tab(m)%rxn_name = 'jcfc115'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CFC-115_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CCl3F (CFC-11) -> Products'
      xsqy_tab(m)%rxn_name = 'jcfcl3'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CFC-11_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => r26
      m = m + 1

      xsqy_tab(m)%equation = 'CCl2F2 (CFC-12) -> Products'
      xsqy_tab(m)%rxn_name = 'jcf2cl2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CFC-12_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => r27
      m = m + 1

      xsqy_tab(m)%equation = 'CH3Br -> Products'
      xsqy_tab(m)%rxn_name = 'jch3br'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CH3Br_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3CCl3 -> Products'
      xsqy_tab(m)%rxn_name = 'jch3ccl3'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CH3CCl3_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => r29
      m = m + 1

      xsqy_tab(m)%equation = 'CH3Cl -> Products'
      xsqy_tab(m)%rxn_name = 'jch3cl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CH3Cl_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => r30
      m = m + 1

      xsqy_tab(m)%equation = 'ClOO -> Products'
      xsqy_tab(m)%rxn_name = 'j_cloo'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/ClOO_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CF3CHCl2 (HCFC-123) -> Products'
      xsqy_tab(m)%rxn_name = 'j_cf3chcl2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      subr(m)%xsqy_sub   => r32
      m = m + 1

      xsqy_tab(m)%equation = 'CF3CHFCl (HCFC-124) -> Products'
      xsqy_tab(m)%rxn_name = 'j_cf3chfcl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      subr(m)%xsqy_sub   => r33
      m = m + 1

      xsqy_tab(m)%equation = 'CH3CFCl2 (HCFC-141b) -> Products'
      xsqy_tab(m)%rxn_name = 'jhcfc141b'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HCFC-141b_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3CF2Cl (HCFC-142b) -> Products'
      xsqy_tab(m)%rxn_name = 'jhcfc142b'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      subr(m)%xsqy_sub   => r35
      m = m + 1

      xsqy_tab(m)%equation = 'CF3CF2CHCl2 (HCFC-225ca) -> Products'
      xsqy_tab(m)%rxn_name = 'j_cf3cf2chcl2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HCFC-225ca_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CF2ClCF2CHFCl (HCFC-225cb) -> Products'
      xsqy_tab(m)%rxn_name = 'j_cf2clcf2chfcl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HCFC-225cb_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CHClF2 (HCFC-22) -> Products'
      xsqy_tab(m)%rxn_name = 'jhcfc22'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HCFC-22_jpl94.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => r38
      m = m + 1

      xsqy_tab(m)%equation = 'HO2 -> OH + O'
      xsqy_tab(m)%rxn_name = 'j_ho2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HO2_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 10
      xsqy_tab(m)%filespec%nread(1) = 15
      subr(m)%xsqy_sub   => r39
      m = m + 1

      xsqy_tab(m)%equation = 'CF2Br2 (Halon-1202) -> Products'
      xsqy_tab(m)%rxn_name = 'j_cf2bf2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Halon-1202_jpl97.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CF2BrCl (Halon-1211) -> Products'
      xsqy_tab(m)%rxn_name = 'jcf2clbr'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Halon-1211_jpl97.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CF3Br (Halon-1301) -> Products'
      xsqy_tab(m)%rxn_name = 'jcf3br'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Halon-1301_jpl97.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CF2BrCF2Br (Halon-2402) -> Products'
      xsqy_tab(m)%rxn_name = 'jh2402'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Halon-2402_jpl97.abs'
      xsqy_tab(m)%filespec%nskip(1) = -1
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'N2O -> N2 + O(1D)'
      xsqy_tab(m)%rxn_name = 'jn2o'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      subr(m)%xsqy_sub   => r44
      m = m + 1

      xsqy_tab(m)%equation   = 'ClONO2 -> Cl + NO3'
      xsqy_tab(m+1)%equation = 'ClONO2 -> ClO + NO2'
      xsqy_tab(m)%rxn_name = 'jclono2_a'
      xsqy_tab(m+1)%rxn_name = 'jclono2_b'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/ClONO2_jpl97.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 119
      subr(m)%xsqy_sub   => r45
      subr(m+1)%xsqy_sub => r45
      m = m + 2

      xsqy_tab(m)%equation   = 'BrONO2 -> BrO + NO2'
      xsqy_tab(m+1)%equation = 'BrONO2 -> Br + NO3'
      xsqy_tab(m)%rxn_name = 'jbrono2_b'
      xsqy_tab(m+1)%rxn_name = 'jbrono2_a'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/BrONO2_jpl03.abs'
      xsqy_tab(m)%filespec%nskip(1) = 13
      xsqy_tab(m)%filespec%nread(1) = 61
      subr(m)%xsqy_sub   => r46
      subr(m+1)%xsqy_sub => r46
      m = m + 2

      xsqy_tab(m)%equation = 'Cl2 -> Cl + Cl'
      xsqy_tab(m)%rxn_name = 'jcl2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      subr(m)%xsqy_sub   => r47
      m = m + 1

      xsqy_tab(m)%equation   = 'HOCH2CHO -> CH2OH + HCO'
      xsqy_tab(m+1)%equation = 'HOCH2CHO -> CH3OH + CO'
      xsqy_tab(m+2)%equation = 'HOCH2CHO -> CH2CHO + OH'
      xsqy_tab(m)%rxn_name = 'j_glyald_a'
      xsqy_tab(m+1)%rxn_name = 'j_glyald_b'
      xsqy_tab(m+2)%rxn_name = 'j_glyald_c'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1:m+2)%channel = (/ 2,3 /)
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CH2OHCHO/glycolaldehyde_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 63
      subr(m)%xsqy_sub   => r101
      subr(m+1)%xsqy_sub => r101
      subr(m+2)%xsqy_sub => r101
      m = m + 3

      xsqy_tab(m)%equation = 'CH3COCOCH3 -> Products'
      xsqy_tab(m)%rxn_name = 'j_biacetyl'
      xsqy_tab(m)%qyld  = .158_rk
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CH3COCOCH3/biacetyl_horowitz.abs'
      xsqy_tab(m)%filespec%nskip(1) = 8
      xsqy_tab(m)%filespec%nread(1) = 287
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3COCH=CH2 -> Products'
      xsqy_tab(m)%rxn_name = 'jmvk'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/MVK_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 146
      subr(m)%xsqy_sub   => r103
      m = m + 1

      xsqy_tab(m)%equation = 'CH2=C(CH3)CHO -> Products'
      xsqy_tab(m)%rxn_name = 'j_macr'
      xsqy_tab(m)%qyld  = .01_rk
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Methacrolein_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 7
      xsqy_tab(m)%filespec%nread(1) = 146
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3COCO(OH) -> Products'
      xsqy_tab(m)%rxn_name = 'j_ch3cocooh'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CH3COCOOH/pyruvic_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 139
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3CH2ONO2 -> CH3CH2O + NO2'
      xsqy_tab(m)%rxn_name = 'j_ch3ch2ono2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/RONO2/RONO2_talukdar.abs'
      xsqy_tab(m)%filespec%nskip(1) = 10
      xsqy_tab(m)%filespec%nread(1) = 63
      subr(m)%xsqy_sub   => r106
      m = m + 1

      xsqy_tab(m)%equation = 'CH3CHONO2CH3 -> CH3CHOCH3 + NO2'
      xsqy_tab(m)%rxn_name = 'j_ch3chono2ch3'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/RONO2/RONO2_talukdar.abs'
      xsqy_tab(m)%filespec%nskip(1) = 10
      xsqy_tab(m)%filespec%nread(1) = 63
      subr(m)%xsqy_sub   => r107
      m = m + 1

      xsqy_tab(m)%equation = 'CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2'
      xsqy_tab(m)%rxn_name = 'j_ch2ohch2ono2'
      xsqy_tab(m)%jndx  = m
      subr(m)%xsqy_sub   => r108
      m = m + 1

      xsqy_tab(m)%equation = 'CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2'
      xsqy_tab(m)%rxn_name = 'j_ch3coch2ono2'
      xsqy_tab(m)%jndx  = m
      subr(m)%xsqy_sub   => r109
      m = m + 1

      xsqy_tab(m)%equation = 'C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2'
      xsqy_tab(m)%rxn_name = 'j_bnit1'
      xsqy_tab(m)%jndx  = m
      subr(m)%xsqy_sub   => r110
      m = m + 1

      xsqy_tab(m)%equation = 'ClOOCl -> Cl + ClOO'
      xsqy_tab(m)%rxn_name = 'j_cloocl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/ClOOCl_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 111
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation   = 'CH2(OH)COCH3 -> CH3CO + CH2(OH)'
      xsqy_tab(m+1)%equation = 'CH2(OH)COCH3 -> CH2(OH)CO + CH3'
      xsqy_tab(m)%rxn_name = 'j_hyac_a'
      xsqy_tab(m+1)%rxn_name = 'j_hyac_b'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Hydroxyacetone_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 96
      subr(m)%xsqy_sub   => r112
      subr(m+1)%xsqy_sub => r112
      m = m + 2

      xsqy_tab(m)%equation = 'HOBr -> OH + Br'
      xsqy_tab(m)%rxn_name = 'jhobr'
      xsqy_tab(m)%jndx  = m
      subr(m)%xsqy_sub   => r113
      m = m + 1 

      xsqy_tab(m)%equation = 'BrO -> Br + O'
      xsqy_tab(m)%rxn_name = 'jbro'
      xsqy_tab(m)%jndx  = m
      subr(m)%xsqy_sub   => r114
      m = m + 1 

      xsqy_tab(m)%equation = 'Br2 -> Br + Br'
      xsqy_tab(m)%rxn_name = 'j_br2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Br2.abs'
      xsqy_tab(m)%filespec%nskip(1) = 6
      xsqy_tab(m)%filespec%nread(1) = 29
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation   = 'NO3-(aq) -> NO2(aq) + O-'
      xsqy_tab(m+1)%equation = 'NO3-(aq) -> NO2-(aq) + O(3P)'
      xsqy_tab(m+2)%equation = 'NO3-(aq) with qy=1'
      xsqy_tab(m)%rxn_name = 'j_no3_aq_a'
      xsqy_tab(m+1)%rxn_name = 'j_no3_aq_b'
      xsqy_tab(m+2)%rxn_name = 'j_no3_aq_c'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1:m+2)%channel = (/ 2,3 /)
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/NO3-_CA03.abs'
      xsqy_tab(m)%filespec%nskip(1) = 7
      xsqy_tab(m)%filespec%nread(1) = 43
      subr(m)%xsqy_sub   => r118
      subr(m+1)%xsqy_sub => r118
      subr(m+2)%xsqy_sub => r118
      m = m + 3

      xsqy_tab(m)%equation = 'CH3COCH2CH3 -> CH3CO + CH2CH3'
      xsqy_tab(m)%rxn_name = 'j_mek'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Martinez.abs'
      xsqy_tab(m)%filespec%nskip(1) = 4
      xsqy_tab(m)%filespec%nread(1) = 96
      subr(m)%xsqy_sub   => r119
      m = m + 1

      xsqy_tab(m)%equation   = 'CH3CH2CO(OONO2) -> CH3CH2CO(OO) + NO2'
      xsqy_tab(m+1)%equation = 'CH3CH2CO(OONO2) -> CH3CH2CO(O) + NO3'
      xsqy_tab(m)%rxn_name = 'j_ppn_a'
      xsqy_tab(m+1)%rxn_name = 'j_ppn_b'
      xsqy_tab(m:m+1)%tpflag  = 1
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/PPN_Harwood.txt'
      xsqy_tab(m)%filespec%nskip(1) = 10
      xsqy_tab(m)%filespec%nread(1) = 66
      subr(m)%xsqy_sub   => r120
      subr(m+1)%xsqy_sub => r120
      m = m + 2

      xsqy_tab(m)%equation = 'HOCH2OOH -> HOCH2O. + OH'
      xsqy_tab(m)%rxn_name = 'j_hoch2ooh'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HOCH2OOH_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 32
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH2=CHCHO -> Products'
      xsqy_tab(m)%rxn_name = 'j_acrol'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Acrolein.txt'
      xsqy_tab(m)%filespec%nskip(1) = 6
      xsqy_tab(m)%filespec%nread(1) = 55
      subr(m)%xsqy_sub   => r122
      m = m + 1

      xsqy_tab(m)%equation = 'CH3CO(OOH) -> Products'
      xsqy_tab(m)%rxn_name = 'j_ch3coooh'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/Peracetic_acid.txt'
      xsqy_tab(m)%filespec%nskip(1) = 6
      xsqy_tab(m)%filespec%nread(1) = 66
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = '(CH3)2NNO -> Products'
      xsqy_tab(m)%rxn_name = 'j_amine'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/dmna.abs'
      xsqy_tab(m)%filespec%nskip(1) = 5
      xsqy_tab(m)%filespec%nread(1) = 132
      xsqy_tab(m)%filespec%xfac(1)  = 1.e-19_rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation   = 'ClO -> Cl + O(1D)'
      xsqy_tab(m+1)%equation = 'ClO -> Cl + O(3P)'
      xsqy_tab(m)%rxn_name = 'j_clo_a'
      xsqy_tab(m+1)%rxn_name = 'jclo'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = 1
      subr(m)%xsqy_sub   => r125
      subr(m+1)%xsqy_sub => r125
      m = m + 2

      xsqy_tab(m)%equation = 'ClNO2 -> Cl + NO2'
      xsqy_tab(m)%rxn_name = 'j_clno2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/ClNO2.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 26
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'BrNO -> Br + NO'
      xsqy_tab(m)%rxn_name = 'j_brno'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/BrNO.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 27
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'BrNO2 -> Br + NO2'
      xsqy_tab(m)%rxn_name = 'j_brno2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/BrNO2.abs'
      xsqy_tab(m)%filespec%nskip(1) = 6
      xsqy_tab(m)%filespec%nread(1) = 54
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation   = 'BrONO -> Br + NO2'
      xsqy_tab(m+1)%equation = 'BrONO -> BrO + NO'
      xsqy_tab(m)%rxn_name = 'jbrono_b'
      xsqy_tab(m+1)%rxn_name = 'jbrono_a'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/BrONO.abs'
      xsqy_tab(m)%filespec%nskip(1) = 8
      xsqy_tab(m)%filespec%nread(1) = 32
      subr(m)%xsqy_sub   => r129
      subr(m+1)%xsqy_sub => r129
      m = m + 2

      xsqy_tab(m)%equation = 'HOCl -> HO + Cl'
      xsqy_tab(m)%rxn_name = 'jhocl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HOCl.abs'
      xsqy_tab(m)%filespec%nskip(1) = 7
      xsqy_tab(m)%filespec%nread(1) = 111
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'NOCl -> NO + Cl'
      xsqy_tab(m)%rxn_name = 'j_nocl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%nfiles = 2
      xsqy_tab(m)%filespec%filename(1:2) = trim(input_data_root)//'/DATAJ1/ABS/NOCl.abs'
      xsqy_tab(m)%filespec%nskip(1:2) = (/ 7,88 /)
      xsqy_tab(m)%filespec%nread(1:2) = (/ 80,61 /)
      subr(m)%xsqy_sub   => r131
      m = m + 1

      xsqy_tab(m)%equation = 'OClO -> Products'
      xsqy_tab(m)%rxn_name = 'joclo'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%nfiles      = 3
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/OClO.abs'
      xsqy_tab(m)%filespec%filename(2) = trim(input_data_root)//'/DATAJ1/ABS/OClO.abs'
      xsqy_tab(m)%filespec%filename(3) = trim(input_data_root)//'/DATAJ1/ABS/OClO.abs'
      xsqy_tab(m)%filespec%nskip(1:3) = (/ 6,1075,2142 /)
      xsqy_tab(m)%filespec%nread(1:3) = (/ 1068,1067,1068 /)
      subr(m)%xsqy_sub   => r132
      m = m + 1

      xsqy_tab(m)%equation = 'BrCl -> Br + Cl'
      xsqy_tab(m)%rxn_name = 'jbrcl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/BrCl.abs'
      xsqy_tab(m)%filespec%nskip(1) = 9
      xsqy_tab(m)%filespec%nread(1) = 81
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3(OONO2) -> CH3(OO) + NO2'
      xsqy_tab(m)%rxn_name = 'j_ch3oono2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CH3OONO2.abs'
      xsqy_tab(m)%filespec%nskip(1) = 9
      xsqy_tab(m)%filespec%nread(1) = 26
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'C(CH3)3(ONO) -> C(CH3)3(O) + NO'
      xsqy_tab(m)%rxn_name = 'j_bnit2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/t-butyl-nitrite.abs'
      xsqy_tab(m)%filespec%nskip(1) = 4
      xsqy_tab(m)%filespec%nread(1) = 96
      xsqy_tab(m)%filespec%xfac(1)  = 1._rk
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'ClONO -> Cl + NO2'
      xsqy_tab(m)%rxn_name = 'j_clono'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/ClONO_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 34
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'HCl -> H + Cl'
      xsqy_tab(m)%rxn_name = 'jhcl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/HCl_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 31
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation   = 'CH2O -> H + HCO' 
      xsqy_tab(m+1)%equation = 'CH2O -> H2 + CO'
      xsqy_tab(m  )%rxn_name = 'jch2o_a'
      xsqy_tab(m+1)%rxn_name = 'jch2o_b'
      xsqy_tab(m)%jndx = m
      xsqy_tab(m+1)%channel = 2
      xsqy_tab(m:m+1)%tpflag = (/ 1,3 /)
      xsqy_tab(m)%filespec%nfiles = 2
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/CH2O/CH2O_jpl11.abs'
      xsqy_tab(m)%filespec%filename(2) = trim(input_data_root)//'/DATAJ1/CH2O/CH2O_jpl11.yld'
      xsqy_tab(m)%filespec%nskip(1:2) = 4
      xsqy_tab(m)%filespec%nread(1:2) = (/ 150,112 /)
      subr(m)%xsqy_sub    => pxCH2O
      subr(m+1)%xsqy_sub  => pxCH2O
      m = m + 2

      xsqy_tab(m)%equation = 'CH3COOH -> CH3 + COOH'
      xsqy_tab(m)%rxn_name = 'j_ch3cooh'
      xsqy_tab(m)%qyld  = .55_rk
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CH3COOH_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 18
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CH3OCl -> CH3O + Cl'
      xsqy_tab(m)%rxn_name = 'j_ch3ocl'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CH3OCl_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 83
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'CHCl3 -> Products'
      xsqy_tab(m)%rxn_name = 'j_chcl3'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/CHCl3_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 39
      subr(m)%xsqy_sub   => r140
      m = m + 1

      xsqy_tab(m)%equation = 'C2H5ONO2 -> C2H5O + NO2'
      xsqy_tab(m)%rxn_name = 'j_c2h5ono2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%tpflag = 1
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/RONO2/C2H5ONO2_iup2006.abs'
      xsqy_tab(m)%filespec%nskip(1) = 4
      xsqy_tab(m)%filespec%nread(1) = 32
      subr(m)%xsqy_sub   => r141
      m = m + 1

      xsqy_tab(m)%equation = 'n-C3H7ONO2 -> C3H7O + NO2'
      xsqy_tab(m)%rxn_name = 'j_nc3h7ono2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/RONO2/nC3H7ONO2_iup2006.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 32
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = '1-C4H9ONO2 -> 1-C4H9O + NO2'
      xsqy_tab(m)%rxn_name = 'j_1c4h9ono2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/RONO2/1C4H9ONO2_iup2006.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 32
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = '2-C4H9ONO2 -> 2-C4H9O + NO2'
      xsqy_tab(m)%rxn_name = 'j_2c4h9ono2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/RONO2/2C4H9ONO2_iup2006.abs'
      xsqy_tab(m)%filespec%nskip(1) = 3
      xsqy_tab(m)%filespec%nread(1) = 15
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'perfluoro 1-iodopropane -> products'
      xsqy_tab(m)%rxn_name = 'j_perfluoro'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/PF-n-iodopropane.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 16
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'I2 -> I + I'
      xsqy_tab(m)%rxn_name = 'j_i2'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/YLD/I2.qy'
      xsqy_tab(m)%filespec%nskip(1) = 4
      xsqy_tab(m)%filespec%nread(1) = 12
      subr(m)%xsqy_sub   => r146
      m = m + 1

      xsqy_tab(m)%equation = 'IO -> I + O'
      xsqy_tab(m)%rxn_name = 'j_io'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/IO_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 133
      subr(m)%xsqy_sub   => no_z_dep
      m = m + 1

      xsqy_tab(m)%equation = 'IOH -> I + OH'
      xsqy_tab(m)%rxn_name = 'j_ioh'
      xsqy_tab(m)%jndx  = m
      xsqy_tab(m)%filespec%filename(1) = trim(input_data_root)//'/DATAJ1/ABS/IOH_jpl11.abs'
      xsqy_tab(m)%filespec%nskip(1) = 2
      xsqy_tab(m)%filespec%nread(1) = 101
      subr(m)%xsqy_sub   => no_z_dep

    end subroutine setup_sub_calls

!-----------------------------------------------------------------------------*
!=  *** ALL the following subroutines have the following arguments
!=  *** except for the routines:
!=      rxn_init, base_read, readit, add_pnts_inter2
!=                                                                           =*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER(len=50) ::, string identifier for each photolysis  (O)=*
!=           reaction defined                                                =*
!-----------------------------------------------------------------------------*

      SUBROUTINE r01(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product of (cross section) x (quantum yield) for the two     =*
!=  O3 photolysis reactions:                                                 =*
!=             (a) O3 + hv -> O2 + O(1D)                                     =*
!=             (b) O3 + hv -> O2 + O(3P)                                     =*
!=  Cross section:  Combined data from WMO 85 Ozone Assessment (use 273K     =*
!=                  value from 175.439-847.5 nm) and data from Molina and    =*
!=                  Molina (use in Hartley and Huggins bans (240.5-350 nm)   =*
!=  Quantum yield:  Choice between                                           =*
!=                   (1) data from Michelsen et al, 1994                     =*
!=                   (2) JPL 87 recommendation                               =*
!=                   (3) JPL 90/92 recommendation (no "tail")                =*
!=                   (4) data from Shetter et al., 1996                      =*
!=                   (5) JPL 97 recommendation                               =*
!=                   (6) JPL 00 recommendation                               =*
!-----------------------------------------------------------------------------*

      use module_xsections, only : o3xs

! input

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      INTEGER, intent(inout) :: j

! local

      INTEGER :: iw
      REAL(rk)    :: xs(nz,nw-1)
      REAL(rk)    :: qy1d(nz)

      errmsg = ' '
      errflg = 0
      xs = xnan
      qy1d = xnan

      if( .not. initialize ) then

! call cross section read/interpolate routine
! cross sections from WMO 1985 Ozone Assessment
! from 175.439 to 847.500 nm. Using value at 273 K.
! Values are over-written in Hartly and Huggins bands, using different
! options depending on value of mopt:

!     mabs = 1 = mostly Reims grp (Malicet, Brion)
!     mabs = 2 = JPL 2006

        CALL o3xs(nz,tlev,nw-1,wl, xs)

!****** quantum yield:
! choose quantum yield recommendation:
!    kjpl87:  JPL recommendation 1987                - JPL 87, 90, 92 do not "tail"
!    kjpl92:  JPL recommendations 1990/92 (identical) - still with no "tail"
!    kjpl97:  JPL recommendation 1997, includes tail, similar to Shetter et al.
!    kmich :  Michelsen et al., 1994
!    kshet :  Shetter et al., 1996
!    kjpl00:  JPL 2000
!    kmats:  Matsumi et al., 2002

! compute cross sections and yields at different wavelengths, altitudes:
        DO iw = 1, nw-1
! quantum yields, Matsumi et al.
          CALL fo3qy2(nz,wc(iw),tlev,qy1d)
          if( xsqy_tab(j)%channel == 2 ) then
            qy1d(1:nz) = (1._rk - qy1d(1:nz))
          endif
          sq(1:nz,iw) = qy1d(1:nz)*xs(1:nz,iw)
        END DO
      endif

      END SUBROUTINE r01

!=============================================================================*

      SUBROUTINE r02(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for NO2            =*
!=  photolysis:                                                              =*
!=         NO2 + hv -> NO + O(3P)                                            =*
!=  Cross section from JPL94 (can also have Davidson et al.)                 =*
!=  Quantum yield from Gardiner, Sperry, and Calvert                         =*
!-----------------------------------------------------------------------------*

      use module_xsections, only : no2xs_jpl06a

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      INTEGER, parameter :: kdata = 200

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata)

! local
      REAL(rk), save :: yg1(kw), ydel(kw)
      REAL(rk) :: yg2(kw)
      REAL(rk) :: qy(nz)
      REAL(rk) :: t(nz)
      REAL(rk) :: no2xs(nz,nw-1)
      INTEGER :: iw, n

      errmsg = ' '
      errflg = 0

 !*************** NO2 photodissociation

      if( initialize ) then
         yg1 = xnan 
         yg2 = xnan
        CALL readit
        ydel(1:nw-1) = yg1(1:nw-1) - yg2(1:nw-1)
      else

! options for NO2 cross section:
! 1 = Davidson et al. (1988), indepedent of T
! 2 = JPL 1994 (same as JPL 1997, JPL 2002)
! 3 = Harder et al.
! 4 = JPL 2006, interpolating between midpoints of bins
! 5 = JPL 2006, bin-to-bin interpolation

!     mabs = 4

        CALL no2xs_jpl06a(nz,tlev,nw,wl, no2xs)

! quantum yields
!     myld = 1   NO2_calvert.yld  (same as JPL2002)
!     myld = 2   NO2_jpl11.yld (same as jpl2006)

!     myld = 2

! from jpl 2011         

        t(1:nz) = .02_rk*(tlev(1:nz) - 298._rk)
        DO iw = 1, nw - 1
          qy(1:nz) = yg1(iw) + ydel(iw)*t(1:nz)
          sq(1:nz,iw) = no2xs(1:nz,iw)*max( qy(1:nz),0._rk )
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 25 ; nsav = 25
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/YLD/NO2_jpl11.yld', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x1,y=y1,y1=y2 )
      xsav(1:n) = x1(1:n)
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/y1(1),0._rk/), errmsg, errflg)
      n = nsav
      x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/y2(1),0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r02

!=============================================================================*

      SUBROUTINE r03(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (absorptioon cross section) x (quantum yield) for    =*
!=  both channels of NO3 photolysis:                                         =*
!=          (a) NO3 + hv -> NO2 + O(3P)                                      =*
!=          (b) NO3 + hv -> NO + O2                                          =*
!=  Cross section combined from Graham and Johnston (<600 nm) and JPL 94     =*
!=  Quantum yield from Madronich (1988)                                      =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      INTEGER, PARAMETER :: kdata=350

      REAL(rk) x(kdata), x1(kdata)
      REAL(rk) y1(kdata)
      real(rk) q1_298(kdata), q1_230(kdata), q1_190(kdata)
      real(rk) q2_298(kdata), q2_230(kdata), q2_190(kdata)
      real(rk) :: sq_wrk(nz)

! local
      real(rk), parameter :: tfac1 = 1._rk/(230._rk - 190._rk)
      real(rk), parameter :: tfac2 = 1._rk/(298._rk - 230._rk)

      REAL(rk) :: xsect
      REAL(rk), save :: yg1(kw)
      real(rk), save :: yg_298(kw,2), yg_230(kw,2), yg_190(kw,2)
      real(rk), save :: delabs(kw,2,2)
      real(rk) :: t(nz)

      INTEGER iw, n, chnl
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
! yields from JPL2011:
          CALL readit
          delabs(1:nw-1,1,1) = yg_230(1:nw-1,1) - yg_190(1:nw-1,1)
          delabs(1:nw-1,2,1) = yg_298(1:nw-1,1) - yg_230(1:nw-1,1)
          delabs(1:nw-1,1,2) = yg_230(1:nw-1,2) - yg_190(1:nw-1,2)
          delabs(1:nw-1,2,2) = yg_298(1:nw-1,2) - yg_230(1:nw-1,2)
          is_initialized = .true.
        endif
      else

! mabs = 3:  JPL11
!     mabs = 3
! myld = 2  from JPL-2011
!     myld = 2

! compute T-dependent quantum yields
        chnl = xsqy_tab(j)%channel
        DO iw = 1, nw-1
          xsect = yg1(iw)
          where(tlev(1:nz) <= 190._rk )
            sq_wrk(1:nz) = yg_190(iw,chnl)*xsect
          elsewhere(tlev(1:nz) > 190._rk .and. tlev(1:nz) <= 230._rk )
            t(1:nz) = tfac1*(tlev(1:nz) - 190._rk)
            sq_wrk(1:nz) = yg_190(iw,chnl) + delabs(iw,1,chnl)*t(1:nz)
          elsewhere(tlev(1:nz) > 230._rk .and. tlev(1:nz) <= 298._rk )
            t(1:nz) = tfac2*(tlev(1:nz) - 230._rk)
            sq_wrk(1:nz) = yg_230(iw,chnl) + delabs(iw,2,chnl)*t(1:nz)
          elsewhere(tlev(1:nz) > 298._rk )
            sq_wrk(1:nz) = yg_298(iw,chnl)
          endwhere
          sq(1:nz,iw) = sq_wrk(1:nz)*xsect
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 289
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/NO3_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=6,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n)*1.E-20_rk
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = 56 ; nsav = 56
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/YLD/NO3_jpl2011.qy', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=5,rd_cnt=n,x=x,y=q1_298, &
                      y1=q1_230,y2=q1_190,y3=q2_298, &
                      y4=q2_230,y5=q2_190 )
      xsav(1:n) = x(1:n)
      q1_298(1:n) = q1_298(1:n)*.001_rk
      q1_230(1:n) = q1_230(1:n)*.001_rk
      q1_190(1:n) = q1_190(1:n)*.001_rk
      q2_298(1:n) = q2_298(1:n)*.001_rk
      q2_230(1:n) = q2_230(1:n)*.001_rk
      q2_190(1:n) = q2_190(1:n)*.001_rk

      CALL add_pnts_inter2(x,q1_298,yg_298,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      n = nsav ; x(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x,q1_230,yg_230,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      n = nsav ; x(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x,q1_190,yg_190,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
     
      n = nsav ; x(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x,q2_298,yg_298(1,2),kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/1._rk,0._rk/), errmsg, errflg)
      n = nsav ; x(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x,q2_230,yg_230(1,2),kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/1._rk,0._rk/), errmsg, errflg)
      n = nsav ; x(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x,q2_190,yg_190(1,2),kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/1._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r03

!=============================================================================*

      SUBROUTINE r04(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yiels) for N2O5 photolysis =*
!=  reactions:                                                               =*
!=       (a) N2O5 + hv -> NO3 + NO + O(3P)                                   =*
!=       (b) N2O5 + hv -> NO3 + NO2                                          =*
!=  Cross section from JPL2011: use tabulated values for 300K, correct for   =*
!=  temperature.
!=  Quantum yield: Analysis of data in JPL94 (->DATAJ1/YLD/N2O5.qy)          =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      INTEGER, PARAMETER :: kdata = 200

      REAL(rk) x1(kdata), x2(kdata)
      REAL(rk) y1(kdata), A(kdata), B(kdata)
      INTEGER :: n1, n2

! local
      INTEGER :: iw
      REAL(rk), save :: yg1(kw), yg2(kw)
      REAL(rk)    :: dum(nz)
      REAL(rk)    :: t(nz)
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else
        if( xsqy_tab(j)%channel == 1 ) then
          DO iw = 1,nw-1
            sq(1:nz,iw) = 0._rk
          ENDDO
        elseif( xsqy_tab(j)%channel == 2 ) then
! temperature dependence only valid for 233 - 295 K.  Extend to 300.
          t(1:nz) = MAX(233._rk,MIN(tlev(1:nz),300._rk))

          DO iw = 1, nw - 1
! Apply temperature correction to 300K values. Do not use A-coefficients 
! because they are inconsistent with the values at 300K.
! quantum yield = 1 for NO2 + NO3, zero for other channels
            dum(1:nz) = 1000._rk*yg2(iw)*(300._rk - t(1:nz))/(300._rk*t(1:nz))
            sq(1:nz,iw) = yg1(iw) * 10._rk**(dum(1:nz))
          ENDDO
        endif
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from jpl2011, at 300 K

      n1 = 103
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/N2O5_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n1,x=x1,y=y1 )
      y1(1:n1) = y1(1:n1) * 1.E-20_rk
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n1, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

! read temperature dependence coefficients:
      n2 = 8
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/N2O5_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=111,rd_cnt=n2,x=x2,y=A,y1=B )

      CALL add_pnts_inter2(x2,B,yg2,kdata,n2, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r04

!=============================================================================*

      SUBROUTINE r06(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yield) for HNO3 photolysis =*
!=        HNO3 + hv -> OH + NO2                                              =*
!=  Cross section: Burkholder et al., 1993                                   =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER n1
      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata)

! local
      real(rk) :: t(nz)
      REAL(rk), save :: yg1(kw), yg2(kw)
      INTEGER i, iw

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else
! quantum yield = 1
! correct for temperature dependence
        t(1:nz) = tlev(1:nz) - 298._rk
        DO iw = 1, nw - 1
          sq(1:nz,iw) = yg1(iw) * exp( yg2(iw)*t(1:nz) )
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
! HNO3 cross section parameters from Burkholder et al. 1993

      integer :: nsav
      real(rk)    :: xsav(kdata)
      real(rk)    :: yends(2)

      n1 =  83 ; nsav = 83
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/HNO3_burk.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=6,rd_cnt=n1,x=y1,y=y2 )

      x1(1:n1) = (/ (184._rk + real(i)*2._rk,i=1,n1) /)
      xsav(1:n1) = x1(1:n1)

      y1(1:n1) = y1(1:n1) * 1.e-20_rk
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n1, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      y2(1:n1) = y2(1:n1) * 1.e-3_rk
      yends(:) = (/ y2(1),y2(n1) /)
      n1 = nsav ; x1(1:n1) = xsav(1:n1)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n1, &
                           nw,wl,xsqy_tab(j)%equation,deltax,yends, errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r06

!=============================================================================*

      SUBROUTINE r08(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yield) for H2O2 photolysis =*
!=         H2O2 + hv -> 2 OH                                                 =*
!=  Cross section:  From JPL97, tabulated values @ 298K for <260nm, T-depend.=*
!=                  parameterization for 260-350nm                           =*
!=  Quantum yield:  Assumed to be unity                                      =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=600

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      real(rk), parameter :: A0 = 6.4761E+04_rk            
      real(rk), parameter :: A1 = -9.2170972E+02_rk        
      real(rk), parameter :: A2 = 4.535649_rk              
      real(rk), parameter :: A3 = -4.4589016E-03_rk        
      real(rk), parameter :: A4 = -4.035101E-05_rk         
      real(rk), parameter :: A5 = 1.6878206E-07_rk
      real(rk), parameter :: A6 = -2.652014E-10_rk
      real(rk), parameter :: A7 = 1.5534675E-13_rk

      real(rk), parameter :: B0 = 6.8123E+03_rk
      real(rk), parameter :: B1 = -5.1351E+01_rk
      real(rk), parameter :: B2 = 1.1522E-01_rk
      real(rk), parameter :: B3 = -3.0493E-05_rk
      real(rk), parameter :: B4 = -1.0924E-07_rk

      INTEGER iw, n
      REAL(rk) lambda
      REAL(rk) sumA, sumB
      REAL(rk) :: t(nz)
      REAL(rk) :: chi(nz)
      REAL(rk), save :: yg(kw)

      errmsg = ' '
      errflg = 0
! cross section from Lin et al. 1978

      if( initialize ) then
        CALL readit
      else
! quantum yield = 1
        t(1:nz) = MIN(MAX(tlev(1:nz),200._rk),400._rk)            
        chi(1:nz) = 1._rk/(1._rk + EXP(-1265._rk/t(1:nz)))
        DO iw = 1, nw - 1
! Parameterization (JPL94)
! Range 260-350 nm; 200-400 K
           IF ((wl(iw) .GE. 260._rk) .AND. (wl(iw) .LT. 350._rk)) THEN
             lambda = wc(iw)
             sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda +  &
                          A4)*lambda +A3)*lambda + A2)*lambda +  &
                          A1)*lambda + A0
             sumB = (((B4*lambda + B3)*lambda + B2)*lambda +  &
                       B1)*lambda + B0

             sq(1:nz,iw) = &
                 (chi(1:nz) * sumA + (1._rk - chi(1:nz))*sumB)*1.E-21_rk
           ELSE
             sq(1:nz,iw) = yg(iw)
           ENDIF
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from JPL94 (identical to JPL97)
! tabulated data up to 260 nm

      integer :: n1

      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/H2O2_jpl94.abs', errmsg=errmsg, errflg=errflg, &
                      rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      
      n1 = 494
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/H2O2_Kahan.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=0,rd_cnt=n1,x=x1(n+1:),y=y1(n+1:) )

      n = n + n1
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r08

!=============================================================================*

      SUBROUTINE r09(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yield) for CHBr3 photolysis=*
!=          CHBr3 + hv -> Products                                           =*
!=  Cross section: Choice of data from Atlas (?Talukdar???) or JPL97         =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=200

      INTEGER n1
      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      REAL(rk), save :: yg(kw)
      real(rk) :: t(nz)

      INTEGER iw

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

! option:

! kopt = 1:  cross section from Elliot Atlas, 1997
! kopt = 2:  cross section from JPL 1997
!     kopt = 2

! quantum yield = 1

        t(1:nz) = 273._rk - tlev(1:nz)
        DO iw = 1, nw - 1
          IF (wc(iw) .GT. 290._rk .AND. wc(iw) .LT. 340._rk ) then
            where( tlev(1:nz) > 210._rk .AND. tlev(1:nz) < 300._rk )
              sq(1:nz,iw) = &
                   EXP( (.06183_rk - .000241_rk*wc(iw))*t(1:nz) &
                             - (2.376_rk + 0.14757_rk*wc(iw)) )
            elsewhere
              sq(1:nz,iw) = yg(iw)
            endwhere
          ELSE
            sq(1:nz,iw) = yg(iw)
          ENDIF
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
! jpl97, with temperature dependence formula,
!w = 290 nm to 340 nm, 
!T = 210K to 300 K
!sigma, cm2 = exp((0.06183-0.000241*w)*(273.-T)-(2.376+0.14757*w))

      n1 = 87
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CHBr3.jpl97', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=6,rd_cnt=n1,x=x1,y=y1 )

      y1(1:n1) = y1(1:n1) * 1.e-20_rk
      CALL add_pnts_inter2(x1,y1,yg,kdata,n1, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/y1(1),0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r09

!=============================================================================*

      SUBROUTINE r11(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: =*
!=      (a)  CH3CHO + hv -> CH3 + HCO                                        =*
!=      (b)  CH3CHO + hv -> CH4 + CO                                         =*
!=      (c)  CH3CHO + hv -> CH3CO + H                                        =*
!=  Cross section:  Choice between                                           =*
!=                   (1) IUPAC 97 data, from Martinez et al.                 =*
!=                   (2) Calvert and Pitts                                   =*
!=                   (3) Martinez et al., Table 1 scanned from paper         =*
!=                   (4) KFA tabulations                                     =*
!=  Quantum yields: Choice between                                           =*
!=                   (1) IUPAC 97, pressure correction using Horowith and    =*
!=                                 Calvert, 1982                             =*
!=                   (2) NCAR data file, from Moortgat, 1986                 =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=150

      INTEGER n
      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata)

! local
      INTEGER :: iw
      INTEGER :: chnl
      REAL(rk)    :: sig
      REAL(rk)    :: qy1_n0, qy1_0, x
      REAL(rk), save :: yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL(rk) :: qy1(nz)
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      chnl = xsqy_tab(j)%channel
      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else
        if( chnl == 1 ) then
!     mabs = 5
!     myld = 1
          DO iw = 1, nw - 1
            sig = yg(iw)
! quantum yields:
! input yields at n0 = 1 atm
            qy1_n0 = yg1(iw)
! Pressure correction for CH3 + CHO channel:
! Assume pressure-dependence only for qy1, not qy2 or qy2.
! Assume total yield 1 at zero pressure
            qy1_0 = 1._rk - (yg2(iw) + yg3(iw))
            
!  compute coefficient:
!  Stern-Volmer:  1/q = 1/q0 + k N  and N0 = 1 atm,
!  then x = K N0 q0 = qy_0/qy_N0 - 1
            if (qy1_n0 > 0._rk) then
              x = qy1_0/qy1_n0 - 1._rk
            else
              x = 0._rk
            endif

            qy1(1:nz) = qy1_n0 * (1._rk + x) / (1._rk + x * airden(1:nz)/2.465E19_rk)
            qy1(1:nz) = MIN( 1._rk,MAX(0._rk,qy1(1:nz)) )
            sq(1:nz,iw) = sig * qy1(1:nz)
          ENDDO
        elseif( chnl == 2 ) then
          sq(1:nw-1,1) = yg(1:nw-1) * yg2(1:nw-1)
        elseif( chnl == 3 ) then
          sq(1:nw-1,1) = yg(1:nw-1) * yg3(1:nw-1)
        endif
      endif

      CONTAINS

      SUBROUTINE readit

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 101
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CH3CHO/CH3CHO_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.e-20_rk

      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

! quantum yields

      n = 12 ; nsav = 12
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CH3CHO/CH3CHO_iup.yld', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n,x=x1,y=y2,y1=y1 )
      xsav(1:n) = x1(1:n)
    
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      n = nsav
      x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      yg3(1:nw-1) = 0._rk

      END SUBROUTINE readit

      END SUBROUTINE r11

!=============================================================================*

      SUBROUTINE r12(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for C2H5CHO        =*
!=  photolysis:                                                              =*
!=         C2H5CHO + hv -> C2H5 + HCO                                        =*
!=                                                                           =*
!=  Cross section:  Choice between                                           =*
!=                   (1) IUPAC 97 data, from Martinez et al.                 =*
!=                   (2) Calvert and Pitts, as tabulated by KFA              =*
!=  Quantum yield:  IUPAC 97 recommendation                                  =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      integer, PARAMETER :: kdata=150

      INTEGER n
      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      REAL(rk), save :: yg(kw), yg1(kw)
      REAL(rk) :: qy1(nz)
      INTEGER iw

      errmsg = ' '
      errflg = 0
      
      if( initialize ) then
        CALL readit
      else

! Absorption:
! 1:  IUPAC-97 data, from Martinez et al.
! 2:  Calvert and Pitts, as tabulated by KFA.

! Quantum yield
! 1:  IUPAC-97 data

!     mabs = 1
!     myld = 1

        DO iw = 1, nw - 1
! quantum yields:
! use Stern-Volmer pressure dependence:
          IF (yg1(iw) .LT. pzero) THEN
            sq(1:nz,iw) = 0._rk
          ELSE
            qy1(1:nz) = 1._rk/(1._rk + (1._rk/yg1(iw) - 1._rk)*airden(1:nz)/2.45e19_rk)
            qy1(1:nz) = MIN(qy1(1:nz),1._rk)
            sq(1:nz,iw) = yg(iw) * qy1(1:nz)
          ENDIF
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      n = 106
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/C2H5CHO/C2H5CHO_iup.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.e-20_rk

      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

! quantum yields

      n = 5
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/C2H5CHO/C2H5CHO_iup.yld', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n,x=x1,y=y1 )

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1._rk-deltax),0._rk,errmsg, errflg)
      CALL addpnt(x1,y1,kdata,n,               0._rk,0._rk,errmsg, errflg)
      CALL addpnt(x1,y1,kdata,n,340._rk,0._rk,errmsg, errflg)
      CALL addpnt(x1,y1,kdata,n,           1.e+38_rk,0._rk,errmsg, errflg)
      CALL inter2(nw,wl,yg1,n,x1,y1,errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r12

!=============================================================================*

      SUBROUTINE r13(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CHOCHO         =*
!=  photolysis:                                                              =*
!=              CHOCHO + hv -> Products                                      =*
!=                                                                           =*
!=  Cross section: Choice between                                            =*
!=                  (1) Plum et al., as tabulated by IUPAC 97                =*
!=                  (2) Plum et al., as tabulated by KFA.                    =*
!=                  (3) Orlando et al.                                       =*
!=                  (4) Horowitz et al., 2001                                =*
!=  Quantum yield: IUPAC 97 recommendation                                   =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=500

      INTEGER n
      REAL(rk) x(kdata), x1(kdata)
      REAL(rk) y1(kdata), y2(kdata), y3(kdata)

! local
      REAL(rk), save :: yg(kw), yg1(kw), yg2(kw), yg3(kw)
      LOGICAL, save :: is_initialized = .false.

!     mabs = 5
!     myld = 2
      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else
        if( xsqy_tab(j)%channel == 1 ) then
          sq(1:nw-1,1) = yg(1:nw-1) * yg1(1:nw-1)
        elseif( xsqy_tab(j)%channel == 2 ) then
          sq(1:nw-1,1) = yg(1:nw-1) * yg2(1:nw-1)
        elseif( xsqy_tab(j)%channel == 3 ) then
          sq(1:nw-1,1) = yg(1:nw-1) * yg3(1:nw-1)
        endif
      endif

      CONTAINS

      SUBROUTINE readit

      integer :: nsav
      real(rk) :: dum(kdata)
      real(rk) :: xsav(kdata)
      real(rk) :: yends(2)

      n = 277
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CHOCHO/glyoxal_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.e-20_rk
      yends(:) = 0._rk
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,yends, errmsg, errflg)

! quantum yields

      n = 40 ; nsav = 40
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CHOCHO/glyoxal_jpl11.qy', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=3,rd_cnt=n,x=x,y=dum,y1=y1,y2=y2,y3=y3 )
      xsav(1:n) = x(1:n)
      yends(1) = y1(1)
      CALL add_pnts_inter2(x,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,yends, errmsg, errflg)
      n = nsav ; x(1:n) = xsav(1:n)
      yends(1) = y2(1)
      CALL add_pnts_inter2(x,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,yends, errmsg, errflg)
      n = nsav ; x(1:n) = xsav(1:n)
      yends(1) = y3(1)
      CALL add_pnts_inter2(x,y3,yg3,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,yends, errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r13

!=============================================================================*

      SUBROUTINE r14(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3COCHO       =*
!=  photolysis:                                                              =*
!=           CH3COCHO + hv -> CH3CO + HCO                                    =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=500

      INTEGER n
      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      REAL(rk), save :: yg(kw)
      REAL(rk) sig
      INTEGER iw
      REAL(rk) phi0, kq

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else
  
!     mabs = 8
!     myld = 5

        DO iw = 1, nw - 1
          sig = yg(iw)
! quantum yields:
! zero pressure yield:
! 1.0 for wc < 380 nm
! 0.0 for wc > 440 nm
! linear in between:
          phi0 = 1._rk - (wc(iw) - 380._rk)/60._rk
          phi0 = MIN(MAX(0._rk,phi0),1._rk)

! Pressure correction: quenching coefficient, torr-1
! in air, Koch and Moortgat:
          kq = 1.36e8_rk * EXP(-8793._rk/wc(iw))
! in N2, Chen et al:
          IF(phi0 .GT. 0._rk) THEN
            IF (wc(iw) .GE. 380._rk .AND. wc(iw) .LE. 440._rk) THEN
              sq(1:nz,iw) = sig * phi0 &
                  / (phi0 + kq * airden(1:nz) * 760._rk/2.456E19_rk)
            ELSE
              sq(1:nz,iw) = sig * phi0
            ENDIF
          ELSE
            sq(1:nz,iw) = 0._rk
          ENDIF
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      n = 294
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CH3COCHO/CH3COCHO_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.e-20_rk
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
         
      END SUBROUTINE readit

      END SUBROUTINE r14

!=============================================================================*

      SUBROUTINE r15(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3COCH3 photolysis=*
!=          CH3COCH3 + hv -> Products                                        =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      integer, PARAMETER :: kdata=150

      INTEGER :: n
      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata), y3(kdata), y4(kdata)

! local
      REAL(rk), save :: yg(kw), yg2(kw), yg3(kw)
      REAL(rk) :: sig(nz)
      REAL(rk) :: T(nz)
      real(rk) :: fac(nz)
      INTEGER iw

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

!     mabs = 4
!     myld = 4

        T(1:nz) = MIN(MAX(tlev(1:nz), 235._rk),298._rk)
        DO iw = 1, nw - 1
          sig(1:nz) = yg(iw) * (1._rk + t(1:nz)*(yg2(iw) + t(1:nz)*yg3(iw)))
          CALL qyacet(nz, wc(iw), tlev, airden, fac)
          sq(1:nz,iw) = sig(1:nz)*min(max(0._rk,fac(1:nz)),1._rk)
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 135 ; nsav = 135
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CH3COCH3/CH3COCH3_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=5,rd_cnt=n,x=x1,y=y1,y1=y2,y2=y3,y3=y4 )
      y1(1:n) = y1(1:n) * 1.e-20_rk
      y2(1:n) = y2(1:n) * 1.e-3_rk
      y3(1:n) = y3(1:n) * 1.e-5_rk
      xsav(1:n) = x1(1:n)

      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y3,yg3,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
         
      END SUBROUTINE readit

      END SUBROUTINE r15

!=============================================================================*

      SUBROUTINE r17(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3ONO2            =*
!=  photolysis:                                                              =*
!=          CH3ONO2 + hv -> CH3O + NO2                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      integer, PARAMETER :: kdata = 100

      INTEGER n
      INTEGER iw
      REAL(rk) :: x1(kdata)
      REAL(rk) :: y1(kdata), y2(kdata)

! local
      REAL(rk), save :: yg(kw), yg1(kw)
      REAL(rk) :: T(nz)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

!     mabs = 9
! quantum yield = 1

        T(1:nz) = tlev(1:nz) - 298._rk
        DO iw = 1, nw - 1
          sq(1:nz,iw) = yg(iw) * exp( yg1(iw) * T(1:nz) )
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 65 ; nsav = 65
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/RONO2/CH3ONO2_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x1,y=y1,y1=y2 )
      y1(1:n) = y1(1:n) * 1.e-20_rk
      y2(1:n) = y2(1:n) * 1.e-3_rk
      xsav(1:n) = x1(1:n)
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r17

!=============================================================================*

      SUBROUTINE r18(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for PAN photolysis:    =*
!=       PAN + hv -> Products                                                =*
!=                                                                           =*
!=  Cross section: from Talukdar et al., 1995                                =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER iw
      INTEGER n
      REAL(rk) :: x1(kdata)
      REAL(rk) :: y1(kdata), y2(kdata)

! local

! quantum yield:
! from JPL 2011 values for >300 nm.
!     real, parameter :: qyNO2 = .7
!     real, parameter :: qyNO3 = .3
      real(rk), parameter :: qyld(2) = (/ .7_rk,.3_rk /)

      INTEGER :: chnl
      REAL(rk), save :: yg(kw), yg2(kw)
      REAL(rk) :: sig(nz), T(nz)
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else
        chnl = xsqy_tab(j)%channel
        T(1:nz) = tlev(1:nz) - 298._rk
        DO iw = 1, nw-1
          sig(1:nz) = yg(iw) * EXP( yg2(iw)*T(1:nz) )
          sq(1:nz,iw) = qyld(chnl) * sig(1:nz)
        ENDDO 
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from 
!      Talukdar et al., 1995, J.Geophys.Res. 100/D7, 14163-14174

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 78 ; nsav = 78
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/RONO2/PAN_talukdar.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=14,rd_cnt=n,x=x1,y=y1,y1=y2 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      y2(1:n) = y2(1:n) * 1.E-3_rk
      xsav(1:n) = x1(1:n)
 
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r18

!=============================================================================*

      SUBROUTINE r20(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CCl4 photolysis:   =*
!=      CCl4 + hv -> Products                                                =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      real(rk), parameter :: b0 = 1.0739_rk
      real(rk), parameter :: b1 = -1.6275e-2_rk
      real(rk), parameter :: b2 = 8.8141e-5_rk
      real(rk), parameter :: b3 = -1.9811e-7_rk
      real(rk), parameter :: b4 = 1.5022e-10_rk

      REAL(rk), save :: yg(kw)
      INTEGER iw, n
      REAL(rk) :: tcoeff
      REAL(rk) :: w1
      REAL(rk) :: temp(nz)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

! mabs = 1:  jpl 1997 recommendation
! mabs = 2:  jpl 2011 recommendation, with T dependence

!     mabs = 2

! compute temperature correction factors:

!** quantum yield assumed to be unity

        temp(1:nz) = min(max(tlev(1:nz),210._rk),300._rk)
        temp(1:nz) = temp(1:nz) - 295._rk
        DO iw = 1, nw-1
! compute temperature correction coefficients:
           tcoeff = 0._rk
           IF(wc(iw) .GT. 194._rk .AND. wc(iw) .LT. 250._rk) THEN 
             w1 = wc(iw)
             tcoeff = b0 + w1*(b1 + w1*(b2 + w1*(b3 + w1*b4)))
           ENDIF
           sq(1:nz,iw) = yg(iw) * 10._rk**(tcoeff*temp(1:nz))
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!** cross sections from JPL97 recommendation (identical to 94 data)

      n = 44
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CCl4_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=5,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
         
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r20

!=============================================================================*

      SUBROUTINE r23(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-113 photolysis:=*
!=          CF2ClCFCl2 + hv -> Products                                      =*
!=  Cross section:  from JPL 97 recommendation, linear interp. between       =*
!=                  values at 210 and 295K                                   =*
!=  Quantum yield:  assumed to be unity                                      =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata)

! local
      real(rk), parameter :: tfac1 = 1._rk/(295._rk - 210._rk)

      REAL(rk), save :: yg2(kw), ydel(kw)
      REAL(rk)       :: yg1(kw)
      REAL(rk) :: t(nz)
      REAL(rk) :: slope(nz)
      INTEGER iw, n

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
        ydel(1:nw-1) = yg1(1:nw-1) - yg2(1:nw-1)
      else

!** quantum yield assumed to be unity

        t(1:nz) = MAX(210._rk,MIN(tlev(1:nz),295._rk))
        slope(1:nz) = (t(1:nz) - 210._rk)*tfac1
        DO iw = 1, nw-1
          sq(1:nz,iw) = yg2(iw) + slope(1:nz)*ydel(iw)
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!** cross sections from JPL97 recommendation (identical to 94 recommendation)

      integer :: nsav
      real(rk)    :: xsav(kdata)

      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CFC-113_jpl94.abs', errmsg=errmsg, errflg=errflg, &
                      rd_cnt=n,x=x1,y=y1,y1=y2 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      y2(1:n) = y2(1:n) * 1.E-20_rk
      xsav(1:n) = x1(1:n)
      nsav = n
      
!* sigma @ 295 K
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

! sigma @ 210 K
      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r23

!=============================================================================*

      SUBROUTINE r24(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-144 photolysis:=*
!=              CF2ClCF2Cl + hv -> Products                                  =*
!=  Cross section: from JPL 97 recommendation, linear interp. between values =*
!=                 at 210 and 295K                                           =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata)

! local
      real(rk), parameter :: tfac1 = 1._rk/(295._rk - 210._rk)

      REAL(rk), save :: yg2(kw), ydel(kw)
      REAL(rk)       :: yg1(kw)
      REAL(rk) :: t(nz)
      REAL(rk) :: slope(nz)
      INTEGER iw, n

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
        ydel(1:nw-1) = yg1(1:nw-1) - yg2(1:nw-1)
      else

!** quantum yield assumed to be unity

        t(1:nz) = MAX(210._rk,MIN(tlev(1:nz),295._rk))
        slope(1:nz) = (t(1:nz) - 210._rk)*tfac1
        DO iw = 1, nw-1
          sq(1:nz,iw) = yg2(iw) + slope(1:nz)*ydel(iw)
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      integer :: nsav
      real(rk)    :: xsav(kdata)

      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CFC-114_jpl94.abs', errmsg=errmsg, errflg=errflg, &
                      rd_cnt=n,x=x1,y=y1,y1=y2 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      y2(1:n) = y2(1:n) * 1.E-20_rk
      xsav(1:n) = x1(1:n)
      nsav = n

!* sigma @ 295 K
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
! sigma @ 210 K
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r24

!=============================================================================*

      SUBROUTINE r26(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-11  photolysis =*
!=          CCl3F + hv -> Products                                           =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      REAL(rk), save :: yg(kw)
      REAL(rk) :: t(nz)
      INTEGER :: iw, n

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

!*** quantum yield assumed to be unity

        t(1:nz) = 1.E-04_rk * (tlev(1:nz) - 298._rk)
        DO iw = 1, nw-1
          sq(1:nz,iw) = yg(iw) * EXP((wc(iw)-184.9_rk) * t(1:nz))
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CFC-11_jpl94.abs', errmsg=errmsg, errflg=errflg,&
                      rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.E-20_rk

!* sigma @ 298 K

      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r26

!=============================================================================*

      SUBROUTINE r27(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-12  photolysis:=*
!=         CCl2F2 + hv -> Products                                           =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      REAL(rk), save :: yg(kw)
      REAL(rk)    :: t(nz)
      INTEGER :: iw, n

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else
!*** quantum yield assumed to be unity
        t(1:nz) = 1.E-04_rk * (tlev(1:nz) - 298._rk) 
        DO iw = 1, nw-1
          sq(1:nz,iw) = yg(iw) * EXP((wc(iw)-184.9_rk) * t(1:nz))
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CFC-12_jpl94.abs', errmsg=errmsg, errflg=errflg, &
                      rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.E-20_rk

!* sigma @ 298 K
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r27

!=============================================================================*

      SUBROUTINE r29(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3CCl3 photolysis =*
!=           CH3CCl3 + hv -> Products                                        =*
!=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
!=                 of data at 210, 250, and 295K                             =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata), y3(kdata)

! local
      real(rk), parameter :: tfac1 = 1._rk/(250._rk - 210._rk)
      real(rk), parameter :: tfac2 = 1._rk/(295._rk - 250._rk)

      REAL(rk), save :: yg2(kw), yg3(kw), ydel1(kw), ydel2(kw)
      REAL(rk)       :: yg1(kw)
      REAL(rk) :: t(nz)
      REAL(rk) :: slope(nz)
      INTEGER iw, n

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
        ydel2(1:nw-1) = yg2(1:nw-1) - yg3(1:nw-1)
        ydel1(1:nw-1) = yg1(1:nw-1) - yg2(1:nw-1)
      else

!*** quantum yield assumed to be unity

        t(1:nz) = MIN(295._rk,MAX(tlev(1:nz),210._rk))
        DO iw = 1, nw-1
          where( t(1:nz) <= 250._rk )
            slope(1:nz) = (t(1:nz) - 210._rk)*tfac1
            sq(1:nz,iw) = yg3(iw) + slope(1:nz)*ydel2(iw)
          elsewhere
            slope(1:nz) = (t(1:nz) - 250._rk)*tfac2
            sq(1:nz,iw) = yg2(iw) + slope(1:nz)*ydel1(iw)
          endwhere
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      integer :: nsav
      real(rk)    :: xsav(kdata)

      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CH3CCl3_jpl94.abs', errmsg=errmsg, errflg=errflg, &
                      rd_cnt=n,x=x1,y=y1,y1=y2,y2=y3 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      y2(1:n) = y2(1:n) * 1.E-20_rk
      y3(1:n) = y3(1:n) * 1.E-20_rk
      xsav(1:n) = x1(1:n)
      nsav = n

!* sigma @ 295 K
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
!* sigma @ 250 K
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
!* sigma @ 210 K
      CALL add_pnts_inter2(x1,y3,yg3,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r29

!=============================================================================*

      SUBROUTINE r30(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3Cl photolysis:  =*
!=            CH3Cl + hv -> Products                                         =*
!=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
!=                 from values at 255, 279, and 296K                         =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata), y3(kdata)

! local
      real(rk), parameter :: tfac1 = 1._rk/(279._rk - 255._rk)
      real(rk), parameter :: tfac2 = 1._rk/(296._rk - 279._rk)

      REAL(rk), save :: yg2(kw), yg3(kw)
      REAL(rk), save :: ydel1(kw), ydel2(kw)
      REAL(rk)       :: yg1(kw)
      REAL(rk) :: t(nz)
      REAL(rk) :: slope(nz)
      INTEGER iw, n

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
        ydel2(1:nw-1) = yg2(1:nw-1) - yg3(1:nw-1)
        ydel1(1:nw-1) = yg1(1:nw-1) - yg2(1:nw-1)
      else

!*** quantum yield assumed to be unity

        t(1:nz) = MAX(255._rk,MIN(tlev(1:nz),296._rk))
        DO iw = 1, nw-1
          where( t(1:nz) <= 279._rk )
            slope(1:nz) = (t(1:nz) - 255._rk)*tfac1
            sq(1:nz,iw) = yg3(iw) + slope(1:nz)*ydel2(iw)
          elsewhere
            slope(1:nz) = (t(1:nz) - 279._rk)*tfac2
            sq(1:nz,iw) = yg2(iw) + slope(1:nz)*ydel1(iw)
          endwhere
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      integer :: nsav
      real(rk)    :: xsav(kdata)

      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CH3Cl_jpl94.abs', errmsg=errmsg, errflg=errflg, &
                      rd_cnt=n,x=x1,y=y1,y1=y2,y2=y3 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      y2(1:n) = y2(1:n) * 1.E-20_rk
      y3(1:n) = y3(1:n) * 1.E-20_rk
      xsav(1:n) = x1(1:n)
      nsav = n

!* sigma @ 296 K
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
!* sigma @ 279 K
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
!* sigma @ 255 K
      CALL add_pnts_inter2(x1,y3,yg3,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r30

!=============================================================================*

      SUBROUTINE r32(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-123 photolysis=*
!=       CF3CHCl2 + hv -> Products                                           =*
!=  Cross section: from Orlando et al., 1991                                 =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
      real(rk), parameter :: LBar = 206.214_rk

      INTEGER i, iw, idum
      INTEGER k
      REAL(rk) lambda
      REAL(rk), save :: TBar
      REAL(rk) :: t(nz)
      REAL(rk) :: sum(nz)
      REAL(rk), save :: coeff(4,3)
      CHARACTER*120 inline

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

!*** quantum yield assumed to be unity

        DO iw = 1, nw-1
          lambda = wc(iw)
! use parameterization only up to 220 nm, as the error bars associated with
! the measurements beyond 220 nm are very large (Orlando, priv.comm.)
          IF (lambda .GE. 190._rk .AND. lambda .LE. 220._rk) THEN
            t(1:nz) = MIN(295._rk,MAX(tlev(1:nz),203._rk)) - TBar
            sum(1:nz) = 0._rk
            DO i = 1, 4
              sum(1:nz) = (coeff(i,1) + t(1:nz)*(coeff(i,2) + t(1:nz)*coeff(i,3))) &
                          * (lambda-LBar)**(i-1) + sum(1:nz)
            ENDDO 
            sq(1:nz,iw) = EXP(sum(1:nz))
          ELSE
            sq(1:nz,iw) = 0._rk
          ENDIF
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross section from Orlando et al., 1991

      OPEN(kin,FILE=trim(input_data_root)//'/DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,'(a120)') inline
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      END SUBROUTINE readit

      END SUBROUTINE r32

!=============================================================================*

      SUBROUTINE r33(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-124 photolysis=*
!=        CF3CHFCl + hv -> Products                                          =*
!=  Cross section: from Orlando et al., 1991                                 =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
      real(rk), parameter :: LBar = 206.214_rk

      INTEGER i, iw, idum
      INTEGER k
      REAL(rk) lambda
      REAL(rk), save :: TBar
      REAL(rk) :: t(nz)
      REAL(rk) :: sum(nz)
      REAL(rk), save :: coeff(4,3)
      CHARACTER*120 inline

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

!*** quantum yield assumed to be unity

        DO iw = 1, nw-1
          lambda = wc(iw)
          IF (lambda .GE. 190._rk .AND. lambda .LE. 230._rk) THEN
            t(1:nz) = MIN(295._rk,MAX(tlev(1:nz),203._rk)) - TBar
            sum(1:nz) = 0._rk
            DO i = 1, 4
              sum(1:nz) = (coeff(i,1) + t(1:nz)*(coeff(i,2) + t(1:nz)*coeff(i,3))) &
                          * (lambda-LBar)**(i-1) + sum(1:nz)
            ENDDO 
            sq(1:nz,iw) = EXP(sum(1:nz))
          ELSE
            sq(1:nz,iw) = 0._rk
          ENDIF
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross section from Orlando et al., 1991

      OPEN(kin,FILE=trim(input_data_root)//'/DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      idum = idum+5
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,'(a120)') inline
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      END SUBROUTINE readit

      END SUBROUTINE r33

!=============================================================================*

      SUBROUTINE r35(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-142b          =*
!=  photolysis:                                                              =*
!=          CH3CF2Cl + hv -> Products                                        =*
!=  Cross section: from Orlando et al., 1991                                 =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
      real(rk), parameter :: LBar = 206.214_rk

      INTEGER i, iw, idum
      INTEGER k
      REAL(rk) lambda
      REAL(rk), save :: Tbar
      REAL(rk) :: t(nz)
      REAL(rk) :: sum(nz)
      REAL(rk), save :: coeff(4,3)
      CHARACTER*80 inline

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

!*** quantum yield assumed to be unity

        DO iw = 1, nw-1
          lambda = wc(iw)
          IF (lambda .GE. 190._rk .AND. lambda .LE. 230._rk) THEN
            t(1:nz) = MIN(295._rk,MAX(tlev(1:nz),203._rk)) - TBar
            sum(1:nz) = 0._rk
            DO i = 1, 4
              sum(1:nz) = (coeff(i,1) + t(1:nz)*(coeff(i,2) + t(1:nz)*coeff(i,3))) &
                          * (lambda-LBar)**(i-1) + sum(1:nz)
            ENDDO 
! offeset exponent by 40 (exp(-40.) = 4.248e-18) to prevent exp. underflow errors
! on some machines.
            sq(1:nz,iw) = 4.248e-18_rk * EXP(sum(1:nz) + 40._rk)
          ELSE
            sq(1:nz,iw) = 0._rk
          ENDIF
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross section from Orlando et al., 1991

      OPEN(kin,FILE=trim(input_data_root)//'/DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      idum = idum+10
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,'(a80)') inline
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      END SUBROUTINE readit

      END SUBROUTINE r35

!=============================================================================*

      SUBROUTINE r38(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-22 photolysis =*
!=          CHClF2 + hv -> Products                                          =*
!=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
!=                 from values at 210, 230, 250, 279, and 295 K              =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

! local
      real(rk), parameter :: tfac1 = 1._rk/(230._rk - 210._rk)
      real(rk), parameter :: tfac2 = 1._rk/(250._rk - 230._rk)
      real(rk), parameter :: tfac3 = 1._rk/(270._rk - 250._rk)
      real(rk), parameter :: tfac4 = 1._rk/(295._rk - 270._rk)

      REAL(rk), save :: yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL(rk)       :: yg1(kw)
      REAL(rk), save :: ydel1(kw), ydel2(kw), ydel3(kw), ydel4(kw)
      REAL(rk) :: t(nz), t1(nz), t2(nz), t3(nz), t4(nz)
      INTEGER iw, n

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
        ydel4(1:nw-1) = yg4(1:nw-1) - yg5(1:nw-1)
        ydel3(1:nw-1) = yg3(1:nw-1) - yg4(1:nw-1)
        ydel2(1:nw-1) = yg2(1:nw-1) - yg3(1:nw-1)
        ydel1(1:nw-1) = yg1(1:nw-1) - yg2(1:nw-1)
      else

!*** quantum yield assumed to be unity

        t(1:nz) = MIN(295._rk,MAX(tlev(1:nz),210._rk))
        t1(1:nz) = (t(1:nz) - 210._rk)*tfac1
        t2(1:nz) = (t(1:nz) - 230._rk)*tfac2
        t3(1:nz) = (t(1:nz) - 250._rk)*tfac3
        t4(1:nz) = (t(1:nz) - 270._rk)*tfac4
        DO iw = 1, nw-1
          where( t(1:nz) <= 230._rk )
            sq(1:nz,iw) = yg5(iw) + t1(1:nz)*ydel4(iw)
          elsewhere( t(1:nz) > 230._rk .and. t(1:nz) <= 250._rk )
            sq(1:nz,iw) = yg4(iw) + t2(1:nz)*ydel3(iw)
          elsewhere( t(1:nz) > 250._rk .and. t(1:nz) <= 270._rk )
            sq(1:nz,iw) = yg3(iw) + t3(1:nz)*ydel2(iw)
          elsewhere
            sq(1:nz,iw) = yg2(iw) + t4(1:nz)*ydel1(iw)
          endwhere
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      integer :: nsav
      real(rk)    :: xsav(kdata)

      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/HCFC-22_jpl94.abs', errmsg=errmsg, errflg=errflg, &
                      rd_cnt=n,x=x1,y=y1,y1=y2,y2=y3,y3=y4,y4=y5 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      y2(1:n) = y2(1:n) * 1.E-20_rk
      y3(1:n) = y3(1:n) * 1.E-20_rk
      y4(1:n) = y4(1:n) * 1.E-20_rk
      y5(1:n) = y5(1:n) * 1.E-20_rk
      nsav = n ; xsav(1:n) = x1(1:n)

!* sigma @ 295 K
      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
!* sigma @ 270 K
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
!* sigma @ 250 K
      CALL add_pnts_inter2(x1,y3,yg3,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
!* sigma @ 230 K
      CALL add_pnts_inter2(x1,y4,yg4,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
!* sigma @ 210 K
      CALL add_pnts_inter2(x1,y5,yg5,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r38

!=============================================================================*

      SUBROUTINE r39(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HO2 photolysis:    =*
!=          HO2 + hv -> OH + O                                               =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed shape based on work by Lee, 1982; normalized      =*
!=                 to unity at 248 nm                                        =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      real(rk), parameter :: tfac1 = 1._rk/(248._rk - 193._rk)
      real(rk), parameter :: xfac1 = 1._rk/15._rk

      REAL(rk), save :: yg(kw)
      REAL(rk), save :: qy(kw)
      INTEGER :: n

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
        WHERE( wc(1:nw-1) >= 248._rk )
          qy(1:nw-1) = 1._rk
        ELSEWHERE
          qy(1:nw-1) = max( (1._rk + (wc(1:nw-1) - 193._rk)*14._rk*tfac1)*xfac1,0._rk )
        ENDWHERE
      else
        sq(1:nw-1,1) = qy(1:nw-1) * yg(1:nw-1)
      endif

      CONTAINS

      SUBROUTINE readit
!*** cross sections from JPL11 recommendation

      n = 15
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/HO2_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=10,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.E-20_rk

      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r39

!=============================================================================*

      SUBROUTINE r44(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for N2O photolysis:    =*
!=              N2O + hv -> N2 + O(1D)                                       =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity, based on Greenblatt and Ravishankara =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
      real(rk), parameter :: A0 = 68.21023_rk                
      real(rk), parameter :: A1 = -4.071805_rk               
      real(rk), parameter :: A2 = 4.301146E-02_rk            
      real(rk), parameter :: A3 = -1.777846E-04_rk           
      real(rk), parameter :: A4 = 2.520672E-07_rk

      real(rk), parameter :: B0 = 123.4014_rk
      real(rk), parameter :: B1 = -2.116255_rk
      real(rk), parameter :: B2 = 1.111572E-02_rk
      real(rk), parameter :: B3 = -1.881058E-05_rk

      INTEGER :: iw
      REAL(rk), save :: a(kw), b(kw)
      REAL(rk) :: lambda
      REAL(rk) :: t(nz)
      REAL(rk) :: bt(nz)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        DO iw = 1, nw-1
          lambda = wc(iw)   
          IF (lambda >= 173._rk .AND. lambda <= 240._rk) THEN
            A(iw) = (((A4*lambda+A3)*lambda+A2)*lambda+A1)*lambda+A0
            B(iw) = (((B3*lambda+B2)*lambda+B1)*lambda+B0)
          ENDIF
        ENDDO
      else

!*** cross sections according to JPL97 recommendation (identical to 94 rec.)
!*** see file DATAJ1/ABS/N2O_jpl94.abs for detail
!*** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
!*** Ravishankara), so quantum yield of O(1D) is assumed to be unity

        t(1:nz) = MAX(194._rk,MIN(tlev(1:nz),320._rk))
        DO iw = 1, nw-1
          lambda = wc(iw)   
          IF (lambda >= 173._rk .AND. lambda <= 240._rk) THEN
            BT(1:nz) = (t(1:nz) - 300._rk)*EXP(B(iw))
            sq(1:nz,iw) = EXP(A(iw)+BT(1:nz))
          ELSE
            sq(1:nz,iw) = 0._rk
          ENDIF
        ENDDO
      endif

      END SUBROUTINE r44

!=============================================================================*

      SUBROUTINE r45(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for ClONO2 photolysis: =*
!=        ClONO2 + hv -> Products                                            =*
!=                                                                           =*
!=  Cross section: JPL 97 recommendation                                     =*
!=  Quantum yield: JPL 97 recommendation                                     =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=150

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata),y2(kdata),y3(kdata)

! local
      REAL(rk) qy1
      REAL(rk) :: xs(nz)
      real(rk) :: t(nz)
      REAL(rk), save :: yg1(kw), yg2(kw), yg3(kw)
      INTEGER iw, chnl
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else

        t(1:nz) = tlev(1:nz) - 296._rk
        chnl = xsqy_tab(j)%channel
        DO iw = 1, nw-1
!** quantum yields (from jpl97, same in jpl2011)
          IF( wc(iw) .LT. 308._rk) THEN
            qy1 = 0.6_rk
          ELSEIF( (wc(iw) .GE. 308) .AND. (wc(iw) .LE. 364._rk) ) THEN
            qy1 = 7.143e-3_rk * wc(iw) - 1.6_rk
          ELSEIF( wc(iw) .GT. 364._rk ) THEN
            qy1 = 1.0_rk
          ENDIF
          IF( chnl == 2 ) then
            qy1 = 1.0_rk - qy1
          ENDIF
! compute T-dependent cross section
          xs(1:nz) = yg1(iw) * (1._rk + t(1:nz) &
                   * (yg2(iw) + t(1:nz)*yg3(iw)))
          sq(1:nz,iw) = qy1 * xs(1:nz)
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
!** cross sections from JPL97 recommendation.  Same in JPL-2011.

      integer :: n
      real(rk)    :: xsav(kz)

      integer, parameter :: nsav = 119

      n = nsav
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/ClONO2_jpl97.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x1,y=y1,y1=y2,y2=y3 )
      xsav(1:n) = x1(1:n)
      y1(1:n)   = y1(1:n) * 1.E-20_rk

      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y3,yg3,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r45

!=============================================================================*

      SUBROUTINE r46(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for BrONO2 photolysis: =*
!=        BrONO2 + hv -> Products                                            =*
!=                                                                           =*
!=  Cross section: JPL 03 recommendation                                     =*
!=  Quantum yield: JPL 03 recommendation                                     =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      REAL(rk), parameter :: qyld(2) = (/ .15_rk,.85_rk /)

      REAL(rk), save :: yg1(kw)
      INTEGER :: n
      INTEGER :: chnl
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else
        chnl = xsqy_tab(j)%channel
        sq(1:nw-1,1) = qyld(chnl) * yg1(1:nw-1)
      endif

      CONTAINS

      SUBROUTINE readit
!** cross sections from JPL03 recommendation

      n = 61
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/BrONO2_jpl03.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=13,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.E-20_rk

      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r46

!=============================================================================*

      SUBROUTINE r47(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for Cl2 photolysis:    =*
!=        Cl2 + hv -> 2 Cl                                                   =*
!=                                                                           =*
!=  Cross section: JPL 97 recommendation                                     =*
!=  Quantum yield: 1     (Calvert and Pitts, 1966)                           =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
      real(rk) :: ex1(nz), ex2(nz)
      real(rk) :: alpha(nz)
      INTEGER iz, iw

      real(rk) :: aa, bb, bb2

      errmsg = ' '
      errflg = 0

      if( .not. initialize ) then

! mabs = 1: Finlayson-Pitts and Pitts
! mabs = 2: JPL2011 formula
      
        DO iz = 1, nz
          aa = 402.7_rk/tlev(iz)
          bb = exp(aa)
          bb2 = bb*bb
          alpha(iz) = (bb2 - 1._rk)/(bb2 + 1._rk)
        ENDDO

!** quantum yield = 1 (Calvert and Pitts, 1966)

        DO iw = 1, nw-1
          ex1(1:nz) = 27.3_rk  * exp(-99.0_rk * alpha(1:nz) * (log(329.5_rk/wc(iw)))**2)
          ex2(1:nz) = .932_rk * exp(-91.5_rk * alpha(1:nz) * (log(406.5_rk/wc(iw)))**2)
          sq(1:nz,iw) = 1.e-20_rk * sqrt(alpha(1:nz)) * (ex1(1:nz) + ex2(1:nz))
        ENDDO
      endif

      END SUBROUTINE r47

!=============================================================================*

      SUBROUTINE r101(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH2(OH)CHO     =*
!=  (glycolaldehye, hydroxy acetaldehyde) photolysis:                        =*
!=           CH2(OH)CHO + hv -> Products                                     =*
!=                                                                           =*
!=  Quantum yield about 50%                                                  =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER :: n
      REAL(rk) x(kdata), y(kdata)

! local
      real(rk), parameter :: qyld(3) = (/ .83_rk, .10_rk, .07_rk /)

      REAL(rk),save :: yg(kw)
      INTEGER :: chnl
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else
        chnl = xsqy_tab(j)%channel
        sq(1:nw-1,1) = yg(1:nw-1) * qyld(chnl)
      endif

      CONTAINS

      SUBROUTINE readit

      n = 63
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CH2OHCHO/glycolaldehyde_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x,y=y )
      y(1:n) = y(1:n) * 1.e-20_rk
         
      CALL add_pnts_inter2(x,y,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r101

!=============================================================================*

      SUBROUTINE r103(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3COCHCH2     =*
!=  Methyl vinyl ketone photolysis:                                          =*
!=           CH3COCH=CH2 + hv -> Products                                     =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=150

      INTEGER n
      REAL(rk) x(kdata), y(kdata)

! local
      REAL(rk), save :: yg(kw)
      REAL(rk) :: qy(nz)
      INTEGER iw

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

! mabs = 1: Schneider and moortgat
! mabs = 2: jpl 2011
!     mabs = 2

! quantum yield from
! Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
! and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
! J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
! depends on pressure and wavelength, set upper limit to 1.0

        DO iw = 1, nw - 1
          qy(1:nz) = exp(-0.055_rk*(wc(iw) - 308._rk)) &
                   / (5.5_rk + 9.2e-19_rk*airden(1:nz))
          qy(1:nz) = min(qy(1:nz), 1._rk)
          sq(1:nz,iw) = yg(iw) * qy(1:nz)
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      n = 146
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/MVK_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x,y=y )
      y(1:n) = y(1:n) * 1.e-20_rk

      CALL add_pnts_inter2(x,y,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r103

!=============================================================================*

      SUBROUTINE r106(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3CH2ONO2     =*
!=  ethyl nitrate       photolysis:                                          =*
!=           CH3CH2ONO2 + hv -> CH3CH2O + NO2                                =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER n1, n2
      REAL(rk) x1(kdata), y1(kdata)
      REAL(rk) x2(kdata), y2(kdata)

! local
      INTEGER iw
      REAL(rk), save :: yg1(kw), yg2(kw)
      real(rk) :: t(nz)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

! quantum yield  = 1

        t(1:nz) = tlev(1:nz) - 298._rk
        DO iw = 1, nw - 1
          sq(1:nz,iw) = yg1(iw)*exp(yg2(iw)*t(1:nz))
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      integer :: n
      real(rk) :: wrk(kdata)

      n = 63
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/RONO2/RONO2_talukdar.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=10,rd_cnt=n,x=x1,y=wrk,y1=wrk, &
                      y2=y1,y3=y2,y4=wrk,y5=wrk )

      x2(1:n) = x1(1:n)

      n1 = count( y1(1:n) > 0._rk )
      if( n1 > 0 ) then
        wrk(1:n1) = pack( y1(1:n),mask=y1(1:n) > 0._rk )
        y1(1:n1)  = wrk(1:n1) * 1.e-20_rk
        wrk(1:n1) = pack( x1(1:n),mask=y1(1:n) > 0._rk )
        x1(1:n1)  = wrk(1:n1)
        CALL add_pnts_inter2(x1,y1,yg1,kdata,n1, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      else
        yg1(:nw) = 0._rk
      endif


      n2 = count( y2(1:n) > 0._rk )
      if( n2 > 0 ) then
        wrk(1:n2) = pack( y2(1:n),mask=y2(1:n) > 0._rk )
        y2(1:n2)  = wrk(1:n2) * 1.e-3_rk
        wrk(1:n2) = pack( x2(1:n),mask=y2(1:n) > 0._rk )
        x2(1:n2)  = wrk(1:n2)
        CALL addpnt(x2,y2,kdata,n2,               0._rk,y2(1),errmsg, errflg)
        CALL addpnt(x2,y2,kdata,n2,           1.e+38_rk,y2(n2),errmsg, errflg)
        CALL inter2(nw,wl,yg2,n2,x2,y2,errmsg, errflg)
      else
        yg2(:nw) = 0._rk
      endif

      END SUBROUTINE readit

      END SUBROUTINE r106

!=============================================================================*

      SUBROUTINE r107(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3CHONO2CH3   =*
!=  isopropyl nitrate   photolysis:                                          =*
!=           CH3CHONO2CH3 + hv -> CH3CHOCH3 + NO2                            =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER n1, n2
      REAL(rk) x1(kdata), y1(kdata)
      REAL(rk) x2(kdata), y2(kdata)

! local
      INTEGER iw
      REAL(rk), save :: yg1(kw), yg2(kw)
      real(rk) :: t(nz)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

! quantum yield  = 1

        t(1:nz) = tlev(1:nz) - 298._rk
        DO iw = 1, nw - 1
          sq(1:nz,iw) = yg1(iw)*exp(yg2(iw)*t(1:nz))
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      integer :: n
      real(rk) :: wrk(kdata)

      n = 63
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/RONO2/RONO2_talukdar.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=10,rd_cnt=n,x=x1,y=wrk, &
                      y1=wrk,y2=wrk,y3=wrk,y4=y1,y5=y2 )

      x2(1:n) = x1(1:n)

      n1 = count( y1(1:n) > 0._rk )
      if( n1 > 0 ) then
        wrk(1:n1) = pack( y1(1:n),mask=y1(1:n) > 0._rk )
        y1(1:n1)  = wrk(1:n1) * 1.e-20_rk
        wrk(1:n1) = pack( x1(1:n),mask=y1(1:n) > 0._rk )
        x1(1:n1)  = wrk(1:n1)
        CALL add_pnts_inter2(x1,y1,yg1,kdata,n1, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      else
        yg1(:nw) = 0._rk
      endif

      n2 = count( y2(1:n) > 0._rk )
      if( n2 > 0 ) then
        wrk(1:n2) = pack( y2(1:n),mask=y2(1:n) > 0._rk )
        y2(1:n2)  = wrk(1:n2) * 1.e-3_rk
        wrk(1:n2) = pack( x2(1:n),mask=y2(1:n) > 0._rk )
        x2(1:n2)  = wrk(1:n2)
        CALL addpnt(x2,y2,kdata,n2,               0._rk,y2(1),errmsg, errflg)
        CALL addpnt(x2,y2,kdata,n2,           1.e+38_rk,y2(n2),errmsg, errflg)
        CALL inter2(nw,wl,yg2,n2,x2,y2,errmsg, errflg)
      else
        yg2(:nw) = 0._rk
      endif

      END SUBROUTINE readit

      END SUBROUTINE r107

!=============================================================================*

      SUBROUTINE r108(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for                =*
!=   nitroxy ethanol CH2(OH)CH2(ONO2) + hv -> CH2(OH)CH2(O.) + NO2           =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
! coefficients from Roberts and Fajer 1989, over 270-306 nm
      real(rk), parameter ::a = -2.359E-3_rk
      real(rk), parameter ::b = 1.2478_rk
      real(rk), parameter ::c = -210.4_rk
      real(rk), save :: xsq(kw)
      
      errmsg = ' '
      errflg = 0

      if( initialize ) then
! quantum yield  = 1
        WHERE( wc(1:nw-1) >= 270._rk .AND. wc(1:nw-1) <= 306._rk )
          xsq(1:nw-1) = EXP(c + wc(1:nw-1)*(b + wc(1:nw-1)*a))
        ELSEWHERE
          xsq(1:nw-1) = 0._rk
        ENDWHERE
      else
        sq(1:nw-1,1) = xsq(1:nw-1)
      endif

      END SUBROUTINE r108

!=============================================================================*

      SUBROUTINE r109(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for                =*
!=   nitroxy acetone CH3COCH2(ONO2) + hv -> CH3COCH2(O.) + NO2               =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
! coefficients from Roberts and Fajer 1989, over 284-335 nm
      real(rk), parameter :: a = -1.365E-3_rk
      real(rk), parameter :: b = 0.7834_rk
      real(rk), parameter :: c = -156.8_rk
      real(rk), save :: xsq(kw)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
! quantum yield  = 1
        WHERE( wc(1:nw-1) >= 284._rk .AND. wc(1:nw-1) <= 335._rk )
          xsq(1:nw-1) = EXP(c + wc(1:nw-1)*(b + wc(1:nw-1)*a))
        ELSEWHERE
          xsq(1:nw-1) = 0._rk
        ENDWHERE
      else
        sq(1:nw-1,1) = xsq(1:nw-1)
      endif

      END SUBROUTINE r109

!=============================================================================*

      SUBROUTINE r110(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for                =*
!=  t-butyl nitrate C(CH3)3(ONO2) + hv -> C(CH3)(O.) + NO2                   =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
! coefficients from Roberts and Fajer 1989, over 270-330 nm
      real(rk), parameter ::a = -0.993E-3_rk
      real(rk), parameter ::b = 0.5307_rk
      real(rk), parameter ::c = -115.5_rk
      real(rk), save :: xsq(kw)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
! quantum yield  = 1
        WHERE( wc(1:nw-1) >= 270._rk .AND. wc(1:nw-1) <= 330._rk )
          xsq(1:nw-1) = EXP(c + wc(1:nw-1)*(b + wc(1:nw-1)*a))
        ELSEWHERE
          xsq(1:nw-1) = 0._rk
        ENDWHERE
      else
        sq(1:nw-1,1) = xsq(1:nw-1)
      endif

      END SUBROUTINE r110

!=============================================================================*

      SUBROUTINE r112(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for hydroxyacetone =*
!=  CH2(OH)COCH3        photolysis:                                          =*
!=           CH2(OH)COCH3  -> CH3CO + CH2OH
!=                         -> CH2(OH)CO + CH3                                =*
!=                                                                           =*
!=  Cross section from Orlando et al. (1999)                                 =*
!=                                                                           =*
!=  Quantum yield assumed 0.325 for each channel (J. Orlando, priv.comm.2003)=*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER :: n
      REAL(rk)    :: x(kdata), y(kdata)

! local
      REAL(rk), parameter :: qy = .325_rk

      REAL(rk) :: yg(kw)
      real(rk), save :: xsq(kw)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
        xsq(1:nw-1) = yg(1:nw-1) * qy
      else
        sq(1:nw-1,1) = xsq(1:nw-1)
      endif

      CONTAINS

      SUBROUTINE readit

      n = 96
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/Hydroxyacetone_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x,y=y )
      y(1:n) = y(1:n) * 1.e-20_rk

      CALL add_pnts_inter2(x,y,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r112

!=============================================================================*

      SUBROUTINE r113(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for HOBr           =*
!=  HOBr -> OH + Br                                                          =*
!=  Cross section from JPL 2003                                              =*
!=  Quantum yield assumed unity as in JPL2003                                =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
      REAL(rk)    :: sig(nw)
      REAL(rk)    :: xfac1(nw)
      real(rk), save :: xsq(kw)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        xsq(1:nw-1) = 0._rk
        WHERE( wc(1:nw-1) >= 250._rk .and. wc(1:nw-1) <= 550._rk )
          xfac1(1:nw-1) = 1._rk/wc(1:nw-1)
          sig(1:nw-1) = 24.77_rk * exp( -109.80_rk*(LOG(284.01_rk*xfac1(1:nw-1)))**2 ) & 
                + 12.22_rk * exp(  -93.63_rk*(LOG(350.57_rk*xfac1(1:nw-1)))**2 ) & 
                + 2.283_rk * exp(- 242.40_rk*(LOG(457.38_rk*xfac1(1:nw-1)))**2 )
          xsq(1:nw-1) = sig(1:nw-1) * 1.e-20_rk
        ENDWHERE
      else
        sq(1:nw-1,1) = xsq(1:nw-1)
      endif

      END SUBROUTINE r113

!=============================================================================*

      SUBROUTINE r114(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for BrO            =*
!=  BrO -> Br + O                                                            =*
!=  Cross section from JPL 2003                                              =*
!=  Quantum yield assumed unity as in JPL2003                                =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! local
      INTEGER :: i, n
      REAL(rk) :: x(20), y(20)
      REAL(rk) :: dum
      REAL(rk) :: yg(kw)
      real(rk), save :: xsq(kw)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAJ1/ABS/BrO.jpl03',STATUS='old')
        DO i = 1, 14
          READ(kin,*)
        ENDDO
        n = 15
        DO i = 1, n
          READ(kin,*) x(i), dum, y(i)
        ENDDO
        CLOSE(kin)

        y(1:n) = y(1:n) * 1.e-20_rk
        n = n + 1
        x(n) = dum
! use bin-to-bin interpolation
        CALL inter4(nw,wl,yg,n,x,y,1, errmsg, errflg)
        xsq(1:nw-1) = yg(1:nw-1)
      else
        sq(1:nw-1,1) = xsq(1:nw-1)
      endif

      END SUBROUTINE r114

!=============================================================================*

      SUBROUTINE r118(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!= NO3-(aq) photolysis for snow simulations                                  =*
!=        a) NO3-(aq) + hv -> NO2 + O-                                       =*
!=        b) NO3-(aq) + hv -> NO2- + O(3P)                                   =*
!=  Cross section:                                                           =*
!=  Burley & Johnston, Geophys. Res. Lett., 19, 1359-1362 (1992)             =*
!=  Chu & Anastasio, J. Phys. Chem. A, 107, 9594-9602 (2003)                 =*
!=  Quantum yield:                                                           =*
!=  Warneck & Wurzinger, J. Phys. Chem., 92, 6278-6283 (1988)                =*
!=  Chu & Anastasio, J. Phys. Chem. A, 107, 9594-9602 (2003)                 =*
!-----------------------------------------------------------------------------*
!= NOTE: user may have to manually add these reactions to the end of the     =*
!= reaction list in file usrinp to include these reactions for a snow run:   =*
!= T 74 NO3-(aq) -> NO2 + O-                                                 =*
!= T 75 NO3-(aq) -> NO2- + O(3P)                                             =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=50

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)    ! y1 = 20'C, y2 = -20'C

! local
      REAL(rk), parameter :: qyld(2:3) = (/ 1.1e-3_rk,1._rk /)
!     REAL, parameter :: qy2 = 1.1e-3
!     REAL, parameter :: qy3 = 1.

      REAL(rk), save :: yg2(kw)
      REAL(rk) :: qy1(nz)
      INTEGER iw, n
      integer :: chnl
      LOGICAL, save :: is_initialized = .false.

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else
        chnl = xsqy_tab(j)%channel
        if( chnl == 1 ) then
          qy1(1:nz) = exp(-2400._rk/tlev(1:nz) + 3.6_rk) ! Chu & Anastasio, 2003
          DO iw = 1, nw-1
            sq(1:nz,iw) = qy1(1:nz)*yg2(iw)
          ENDDO
        else
          sq(1:nw-1,1) = qyld(chnl)*yg2(1:nw-1)
        endif
      endif

      CONTAINS

      SUBROUTINE readit
!** NO3-(aq) cross sections from Chu and Anastasio 2003:
! convert from molar abs log10 to cm2 per molec

      real(rk) :: wrk(kdata)

      n = 43
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/NO3-_CA03.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=7,rd_cnt=n,x=x1,y=y1,y1=wrk, &
                      y2=wrk,y3=wrk,y4=wrk )
      y1(1:n) = y1(1:n) * 3.82e-21_rk
      CALL add_pnts_inter2(x1,y1,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r118

!=============================================================================*

      SUBROUTINE r119(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for                =*
!=    methylethylketone                                                      =*
!=  CH3COCH2CH3 photolysis:                                                  =*
!=           CH3COCH2CH3  -> CH3CO + CH2CH3                                  =*
!=                                                                           =*
!=  Cross section from Martinez et al. (1992)                                =*
!=                                                                           =*
!=  Quantum yield assumed 0.325 for each channel (J. Orlando, priv.comm.2003)=*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER n
      REAL(rk) x(kdata), y(kdata)

! local
      REAL(rk), save :: yg(kw)
      REAL(rk) :: ptorr(nz)
      REAL(rk) :: qy(nz)
      INTEGER iw

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

! Quantum Yields from 
! Raber, W.H. (1992) PhD Thesis, Johannes Gutenberg-Universitaet, Mainz, Germany.
! other channels assumed negligible (less than 10%).
! Total quantum yield  = 0.38 at 760 Torr.
! Stern-Volmer form given:  1/phi = 0.96 + 2.22e-3*P(torr)
!     compute local pressure in torr

        ptorr(1:nz) = 760._rk*airden(1:nz)/2.69e19_rk
        qy(1:nz)    = min( 1._rk/(0.96_rk + 2.22E-3_rk*ptorr(1:nz)),1._rk )
        DO iw = 1, nw-1
          sq(1:nz,iw) = yg(iw) * qy(1:nz)
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      real(rk) :: wrk(kdata)
      n = 96
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/Martinez.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n,x=x,y=wrk,y1=y, &
                      y2=wrk,y3=wrk )
      y(1:n) = y(1:n) * 1.e-20_rk

      CALL add_pnts_inter2(x,y,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r119

!=============================================================================*

      SUBROUTINE r120(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for PPN photolysis:    =*
!=       PPN + hv -> Products                                                =*
!=                                                                           =*
!=  Cross section: from JPL 2006 (originally from Harwood et al. 2003)       =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER :: iw
      INTEGER :: n
      REAL(rk)    :: x1(kdata)
      REAL(rk)    :: y1(kdata), y2(kdata)

! local
      real(rk), parameter :: qyld(2) = (/ 0.61_rk,0.39_rk /)

      INTEGER :: chnl
      REAL(rk), save :: yg(kw), yg2(kw)
      real(rk) :: t(nz)
      REAL(rk) :: sig(nz)
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else
    
        chnl = xsqy_tab(j)%channel
        t(1:nz) = tlev(1:nz) - 298._rk
        DO iw = 1, nw-1
          sig(1:nz) = yg(iw) * EXP(yg2(iw)*t(1:nz))
          sq(1:nz,iw) = qyld(chnl) * sig(1:nz)
        ENDDO 
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from JPL 2011 (originally from Harwood et al. 2003)

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 66 ; nsav = 66
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/PPN_Harwood.txt', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=10,rd_cnt=n,x=x1,y=y1,y1=y2 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      y2(1:n) = y2(1:n) * 1E-3_rk
      xsav(1:n) = x1(1:n)
 
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r120

!=============================================================================*

      SUBROUTINE r122(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH2=CHCHO          =*
!=  (acrolein) photolysis:                                                   =*
!=       CH2=CHCHO + hv -> Products                                          =*
!=                                                                           =*
!=  Cross section: from JPL 2006 (originally from Magneron et al.            =*
!=  Quantum yield: P-dependent, JPL 2006 orig. from Gardner et al.           =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=100

      INTEGER iw
      INTEGER n
      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      REAL(rk), save :: yg(kw)
      real(rk) :: qy(nz), qym1(nz)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else

! quantum yields are pressure dependent between air number densities
! of 8e17 and 2.6e19, Gardner et al.:
        DO iw = 1, nw-1
          where( airden(1:nz) > 2.6e19_rk )
            qy(1:nz) = 0.004_rk
          elsewhere( airden(1:nz) > 8.e17_rk .and. airden(1:nz) <= 2.6e19_rk )
            qym1(1:nz) = 0.086_rk + 1.613e-17_rk * airden(1:nz)
            qy(1:nz)   = 0.004_rk + 1._rk/qym1(1:nz)
          elsewhere( airden(1:nz) <= 8.e17_rk )
            qym1(1:nz) = 0.086_rk + 1.613e-17_rk * 8.e17_rk
            qy(1:nz)   = 0.004_rk + 1._rk/qym1(1:nz)
          endwhere
          sq(1:nz,iw) = qy(1:nz) * yg(iw)
        ENDDO 
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from JPL 2006 (originally from Magneron et al.)

      n = 55
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/Acrolein.txt', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=6,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
 
      CALL add_pnts_inter2(x1,y1,yg,kdata,n,nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r122

!=============================================================================*

      SUBROUTINE r125(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for ClO photolysis     =*
!=       ClO + hv -> Cl + O                                                  =*
!=                                                                           =*
!=  Cross section: from Maric and Burrows 1999                               =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      INTEGER, intent(inout) :: j
! data arrays
      integer, PARAMETER :: kdata=500

      INTEGER iw
      INTEGER i
      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
      REAL(rk) :: yg(kw)
      REAL(rk) qy1, qy2

      real(rk), save :: tmp(12)
      real(rk), save :: ygt(kw,12)
      real(rk) x(kdata), y(kdata,12)
      real(rk) tx, xdum
      integer m, nn, ii
      real(rk) yy
      INTEGER m1, m2
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          tmp(1)    = 180._rk
          tmp(2:12) = (/ (190._rk + 10._rk*real(m-1),m=2,12) /)
          is_initialized = .true.
        endif
      else

        DO i = 1, nz
          tx = tlev(i)
! locate temperature indices for interpolation:
          m1 = 1 + INT(.1_rk*(tx - 190._rk))
          m1 = MIN(MAX(1 ,m1),11)
          m2 = m1 + 1
          DO iw = 1, nw-1
            yy = ygt(iw,m1) + (ygt(iw,m2) - ygt(iw,m1)) &
                 * (tx - tmp(m1))/(tmp(m2) - tmp(m1))
! threshold for O(1D) productionis 263.4 nm:
            if(wc(iw) .lt. 263.4_rk) then
               qy1 = 1._rk
            else
               qy1 = 0._rk
            endif
            qy2 = 1._rk - qy1
            if( xsqy_tab(j)%channel == 1 ) then
              sq(i,iw) = qy1 * yy
            elseif( xsqy_tab(j)%channel == 2 ) then
              sq(i,iw) = qy2 * yy
            endif
          ENDDO
        ENDDO 
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from 
! Maric D. and J.P. Burrows, J. Quantitative Spectroscopy and 
! Radiative Transfer 62, 345-369, 1999.  Data was downloaded from 
! their web site on 15 September 2009.

      integer :: nsav
      real(rk)    :: xsav(kdata)

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAJ1/ABS/ClO_spectrum.prn',STATUS='OLD')
      DO i = 1, 2
         READ(kin,*)
      ENDDO
      nn = 453 ; nsav = 453
      DO ii = 1, nn
         i = nn - ii + 1
         READ(kin,*) xdum, x(i), xdum, (y(i,m), m = 1, 12)
      ENDDO
      CLOSE(kin)

      xsav(1:nn) = x(1:nn)
      DO m = 1, 12
         nn = nsav
         x1(1:nn) = xsav(1:nn)
         y1(1:nn) = y(1:nn,m)
         CALL add_pnts_inter2(x1,y1,yg,kdata,nn, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
         ygt(1:nw-1,m) = yg(1:nw-1)
      ENDDO

      END SUBROUTINE readit

      END SUBROUTINE r125

!=============================================================================*

      SUBROUTINE r129(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for bromine nitrite    =*
!=       BrONO -> Br + NO2                                                   =*
!=       BrONO -> BrO + NO                                                   =*
!=                                                                           =*
!=  Cross section: from IUPAC (vol.3)                                        =*
!=  Quantum yield: Assumed to be 0.5 for each channel                        =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=50

      INTEGER :: n
      INTEGER :: chnl
      REAL(rk)    :: x1(kdata)
      REAL(rk)    :: y1(kdata)

! local
      real(rk), parameter :: qyld(2) = 0.5_rk

      REAL(rk) :: yg(kw)
      real(rk), save :: xsq(kw)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
        chnl = xsqy_tab(j)%channel
        xsq(1:nw-1) = qyld(chnl) * yg(1:nw-1)
      else
        sq(1:nw-1,1) = xsq(1:nw-1)
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from IUPAC (vol III) 2007

      n = 32
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/BrONO.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=8,rd_cnt=n,x=x1,y=y1 )
 
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r129

!******************************************************************

      SUBROUTINE r131(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for 
!=       NOCl -> NO + Cl                                                     =*
!=  Cross section: from IUPAC (vol.3)                                        =*
!=  Quantum yield: Assumed to be 1                                           =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=150

      INTEGER iw
      INTEGER n, ii
      REAL(rk) x1(kdata), y1(kdata)
      REAL(rk) y223(kdata),y243(kdata),y263(kdata),y298(kdata), &
           y323(kdata), y343(kdata)

! local
      REAL(rk), save :: yg223(kw),yg243(kw),yg263(kw), &
                    yg298(kw),yg323(kw), yg343(kw)
      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else
! quantum yields assumed unity
        DO iw = 1, nw-1
          where( tlev(1:nz) .le. 223._rk )
            sq(1:nz,iw) = yg223(iw)
          elsewhere (tlev(1:nz) .gt. 223._rk .and. tlev(1:nz) .le. 243._rk )
            sq(1:nz,iw) = yg223(iw) &
                   + (yg243(iw) - yg223(iw))*(tlev(1:nz) - 223._rk)*.05_rk
          elsewhere (tlev(1:nz) .gt. 243._rk .and. tlev(1:nz) .le. 263._rk )
            sq(1:nz,iw) = yg243(iw) &
                   + (yg263(iw) - yg243(iw))*(tlev(1:nz) - 243._rk)*.05_rk
          elsewhere (tlev(1:nz) .gt. 263._rk .and. tlev(1:nz) .le. 298._rk )
            sq(1:nz,iw) = yg263(iw) &
                   + (yg298(iw) - yg263(iw))*(tlev(1:nz) - 263._rk)/35._rk
          elsewhere (tlev(1:nz) .gt. 298._rk .and. tlev(1:nz) .le. 323._rk )
            sq(1:nz,iw) = yg298(iw) &
                   + (yg323(iw) - yg298(iw))*(tlev(1:nz) - 298._rk)*.04_rk
          elsewhere (tlev(1:nz) .gt. 323._rk .and. tlev(1:nz) .le. 343._rk )
            sq(1:nz,iw) = yg323(iw) &
                   + (yg343(iw) - yg323(iw))*(tlev(1:nz) - 323._rk)*.05_rk
          elsewhere (tlev(1:nz) .gt. 343._rk )
            sq(1:nz,iw) = 0._rk
          endwhere
        ENDDO 
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from IUPAC (vol III) 2007

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 80
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/NOCl.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=7,rd_cnt=n,x=x1,y=y1 )
      y223(1:n) = y1(1:n)
      y243(1:n) = y1(1:n)
      y263(1:n) = y1(1:n)
      y298(1:n) = y1(1:n)
      y323(1:n) = y1(1:n)
      y343(1:n) = y1(1:n)
      ii = 61
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/NOCl.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=88,rd_cnt=ii,x=x1(n+1:),y=y223(n+1:), &
                      y1=y243(n+1:),y2=y263(n+1:),y3=y298(n+1:), &
                      y4=y323(n+1:),y5=y343(n+1:) )
      
      n = n + ii
      nsav = n ; xsav(1:n) = x1(1:n)

      CALL add_pnts_inter2(x1,y223,yg223,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y243,yg243,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y263,yg263,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y298,yg298,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y323,yg323,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y343,yg343,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r131

!******************************************************************

      SUBROUTINE r132(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for 
!=       OClO -> Products                                                    =*
!=  Cross section: from Wahner et al., J. Phys. Chem. 91, 2734, 1987         =*
!=  Quantum yield: Assumed to be 1                                           =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=2000

      INTEGER iw
      INTEGER i
      integer n204, n296, n378
      REAL(rk) x204(kdata),x296(kdata),x378(kdata)
      REAL(rk) y204(kdata),y296(kdata),y378(kdata)

! local
      REAL(rk), save :: yg204(kw),yg296(kw),yg378(kw)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else
! quantum yields assumed unity
        DO iw = 1, nw-1
          where(tlev(1:nz) .le. 204._rk )
            sq(1:nz,iw) = yg204(iw)
          elsewhere (tlev(1:nz) .gt. 204._rk .and. tlev(1:nz) .le. 296._rk )
            sq(1:nz,iw) = yg204(iw) &
                + (yg296(iw) - yg204(iw))*(tlev(1:nz) - 204._rk)/92._rk
          elsewhere (tlev(1:nz) .gt. 296._rk .and. tlev(1:nz) .le. 378._rk )
            sq(1:nz,iw) = yg296(iw) &
                + (yg378(iw) - yg296(iw))*(tlev(1:nz) - 296._rk)/82._rk
          elsewhere (tlev(1:nz) .gt. 378._rk )
            sq(1:nz,iw) = yg378(iw)  
          endwhere
        ENDDO 
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from 
!A. Wahner, G.S. tyndall, A.R. Ravishankara, J. Phys. Chem., 91, 2734, (1987).
!Supplementary Data, as quoted at:
!http://www.atmosphere.mpg.de/enid/26b4b5172008b02407b2e47f08de2fa1,0/Spectra/Introduction_1rr.html

      OPEN(UNIT=kin,FILE=trim(input_data_root)//'/DATAJ1/ABS/OClO.abs',STATUS='OLD')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n204 = 1074-6
      DO i = 1, n204
         READ(kin,*) x204(i), y204(i)
      ENDDO

      READ(kin,*)
      n296 = 1067
      do i = 1, n296
         read(kin,*) x296(i), y296(i)
      enddo

      read(kin,*)
      n378 = 1068
      do i = 1, n378
         read(kin,*) x378(i), y378(i)
      enddo

      CLOSE(kin)
      
      CALL add_pnts_inter2(x204,y204,yg204,kdata,n204, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      CALL add_pnts_inter2(x296,y296,yg296,kdata,n296, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      CALL add_pnts_inter2(x378,y378,yg378,kdata,n378, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r132

!******************************************************************

      SUBROUTINE pxCH2O(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  JPL 2011 recommendation.                                                 =*
!=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
!=        (a) CH2O + hv -> H + HCO                                           =*
!=        (b) CH2O + hv -> H2 + CO                                           =*
!=  written by s. madronich march 2013
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      integer, PARAMETER :: kdata=200

! data arrays
      INTEGER iw
      INTEGER n
      REAL(rk) x1(kdata)
      REAL(rk) y298(kdata), tcoef(kdata)
      REAL(rk) qr(kdata), qm(kdata)

! local
      REAL(rk) ak300
      real(rk) qyr300, qym300
      REAL(rk), save :: yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL(rk) :: t(nz), t1(nz)
      REAL(rk) :: sig(nz)
      REAL(rk) :: qymt(nz)
      REAL(rk) :: akt(nz)
      LOGICAL, save :: is_initialized = .false.

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        if( .not. is_initialized ) then
          CALL readit
          is_initialized = .true.
        endif
      else

        t(1:nz)  = tlev(1:nz) - 298._rk
        t1(1:nz) = (300._rk - tlev(1:nz))/80._rk
        DO iw = 1, nw - 1
! correct cross section for temperature dependence:
          sig(1:nz) = yg1(iw) + yg2(iw) * t(1:nz)
! assign room temperature quantum yields for radical and molecular channels
          qyr300 = yg3(iw)
          qym300 = yg4(iw)
! between 330 ande 360 nm, molecular channel is pressure and temperature dependent.
          IF (wc(iw) .ge. 330._rk .and. wc(iw) .lt. 360._rk .and. qym300 .gt. 0._rk) then
            ak300 = (1._rk - (qym300+qyr300))/(qym300*(1._rk - qyr300))
            ak300 = ak300/2.45e19_rk
            akt(1:nz) = ak300 * (1._rk + 0.05_rk * (wc(iw) - 329._rk) * t1(1:nz))
            qymt(1:nz) = 1._rk/(1._rk/(1._rk-qyr300) + akt(1:nz)*airden(1:nz))
          ELSE
            qymt(1:nz) = qym300
          ENDIF
          if( xsqy_tab(j)%channel == 1 ) then
            sq(1:nz,iw) = sig(1:nz) * qyr300
          elseif( xsqy_tab(j)%channel == 2 ) then
            sq(1:nz,iw) = sig(1:nz) * qymt(1:nz)
          endif
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
! read JPL2011 cross section data:

      integer :: nsav
      real(rk)    :: xsav(kdata)

      n = 150 ; nsav = 150
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CH2O/CH2O_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n,x=x1,y=y298, &
                      y1=tcoef )
      y298(1:n)  = y298(1:n) * 1.e-20_rk
      tcoef(1:n) = tcoef(1:n) * 1.e-24_rk
      xsav(1:n) = x1(1:n)
      
!     terminate endpoints and interpolate to working grid
      CALL add_pnts_inter2(x1,y298,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      
      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,tcoef,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

! quantum yields: Read, terminate, interpolate:

      n = 112 ; nsav = 112
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/CH2O/CH2O_jpl11.yld', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n,x=x1,y=qr,y1=qm )
      xsav(1:n) = x1(1:n)

      CALL add_pnts_inter2(x1,qr,yg3,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/qr(1),0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,qm,yg4,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/qm(1),0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE pxCH2O

!=============================================================================*

      SUBROUTINE r140(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CHCl3 photolysis:  =*
!=      CHCL3 + hv -> Products                                               =*
!=  Cross section: from JPL 2011 recommendation                              =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=50

      REAL(rk) x1(kdata)
      REAL(rk) y1(kdata)

! local
! temperature correction factors:
      real(rk), parameter :: b0 = 3.7973_rk
      real(rk), parameter :: b1 = -7.0913e-2_rk
      real(rk), parameter :: b2 = 4.9397e-4_rk
      real(rk), parameter :: b3 = -1.5226e-6_rk
      real(rk), parameter :: b4 = 1.7555e-9_rk

      INTEGER :: iw, n
      REAL(rk), save :: yg(kw)
      REAL(rk)    :: tcoeff
      REAL(rk)    :: w1
      REAL(rk)    :: temp(nz)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else
      
!** quantum yield assumed to be unity
        temp(1:nz) = min(max(tlev(1:nz),210._rk),300._rk) - 295._rk
        DO iw = 1, nw-1
! compute temperature correction coefficients:
          tcoeff = 0._rk
          w1 = wc(iw)
          IF(w1 > 190._rk .AND. w1 < 240._rk) THEN 
            tcoeff = b0 + w1*(b1 + w1*(b2 + w1*(b3 + w1*b4)))
          ENDIF
          sq(1:nz,iw) = yg(iw) * 10._rk**(tcoeff*temp(1:nz))
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit

      n = 39
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/CHCl3_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=3,rd_cnt=n,x=x1,y=y1 )
      y1(1:n) = y1(1:n) * 1.E-20_rk
      
      CALL add_pnts_inter2(x1,y1,yg,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r140

!=============================================================================*

      SUBROUTINE r141(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for C2H5ONO2           =*
!=  photolysis:                                                              =*
!=          C2H5ONO2 + hv -> C2H5O + NO2                                     =*
!=                                                                           =*
!=  Cross section:  IUPAC 2006 (Atkinson et al., ACP, 6, 3625-4055, 2006)    =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata = 50

      INTEGER :: iw
      REAL(rk)    :: x1(kdata)
      REAL(rk)    :: y1(kdata), y2(kdata)

! local
      REAL(rk), save :: yg1(kw), yg2(kw)
      real(rk) :: t(nz)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else
! quantum yield = 1
        t(1:nz) = tlev(1:nz) - 298._rk
        DO iw = 1, nw - 1
          sq(1:nz,iw) = yg1(iw) * exp(yg2(iw) * t(1:nz))
        ENDDO
      endif

      CONTAINS

      SUBROUTINE readit
! mabs: absorption cross section options: 1:  IUPAC 2006

      integer :: n, nsav
      real(rk)    :: xsav(kdata)

      n = 32 ; nsav = 32
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/RONO2/C2H5ONO2_iup2006.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n,x=x1,y=y1,y1=y2 )
      y1(1:n) = y1(1:n) * 1.e-20_rk
      y2(1:n) = y2(1:n) * 1.e-3_rk
      xsav(1:n) = x1(1:n)

      CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      n = nsav ; x1(1:n) = xsav(1:n)
      CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r141

      SUBROUTINE r146(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for                    =*
!=     molecular Iodine, I2                                                  =*
!=  cross section from JPL2011                                               =*
!=  Quantum yield: wave-dep, from Brewer and Tellinhuisen, 1972              =*
!=  Quantum yield for Unimolecular Dissociation of I2 in Visible Absorption  =*
!=  J. Chem. Phys. 56, 3929-3937, 1972.
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

! data arrays
      integer, PARAMETER :: kdata=200

      INTEGER :: n
      REAL(rk)    :: x(kdata), y(kdata)

! local
      REAL(rk), save    :: yg1(kw), yg2(kw)

      errmsg = ' '
      errflg = 0

      if( initialize ) then
        CALL readit
      else  
        sq(1:nw-1,1) = yg1(1:nw-1) * yg2(1:nw-1)
      endif

      CONTAINS

      SUBROUTINE readit
! cross section from JPL2011

      n = 104
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/ABS/I2_jpl11.abs', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=2,rd_cnt=n,x=x,y=y )
      y(1:n) = y(1:n) * 1.e-20_rk
      
      CALL add_pnts_inter2(x,y,yg1,kdata,n, &
                           nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

! quantum yields 

      n = 12
      CALL base_read( filespec=trim(input_data_root)//'/DATAJ1/YLD/I2.qy', errmsg=errmsg, errflg=errflg, &
                      skip_cnt=4,rd_cnt=n,x=x,y=y )
      
      CALL add_pnts_inter2(x,y,yg2,kdata,n,nw,wl,xsqy_tab(j)%equation,deltax,(/1._rk,0._rk/), errmsg, errflg)

      END SUBROUTINE readit

      END SUBROUTINE r146

      subroutine XSQY_O3(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide the product of (cross section) x (quantum yield) for the two      !
!   o3 photolysis reactions:                                                  !
!              (a) O3 + hv -> O2 + O(1D)                                      !
!              (b) O3 + hv -> O2 + O(3P)                                      !
!   cross sections:                                                           !
!               120nm - 185 nm = Ackerman, 1971                               !
!               185nm - 827 nm = JPL06 (293-298K)                             !
!               196nm - 342 nm = JPL06 (218 K)                                !
!   quantum yield:  JPL 06 recommendation                                     !
!-----------------------------------------------------------------------------!
!   parameters:                                                               !
!   nw     - integer, number of specified intervals + 1 in working        (i) !
!            wavelength grid                                                  !
!   wl     - real, vector of lower limits of wavelength intervals in      (i) !
!            working wavelength grid                                          !
!   wc     - real, vector of center points of wavelength intervals in     (i) !
!            working wavelength grid                                          !
!   nz     - integer, number of altitude levels in working altitude grid  (i) !
!   tlev   - real, temperature (k) at each specified altitude level       (i) !
!   airlev - real, air density (molec/cc) at each altitude level          (i) !
!   j      - integer, counter for number of weighting functions defined  (io) !
!   sq     - real, cross section x quantum yield (cm^2) for each          (o) !
!            photolysis reaction defined, at each defined wavelength and      !
!            at each defined altitude level                                   !
!   jlabel - character*60, string identifier for each photolysis reaction (o) !
!            defined                                                          !
!-----------------------------------------------------------------------------!
!   edit history:                                                             !
!   02/06/08  JPL-06 Doug Kinnison                                            !
!   01/02  Atkinson, 1971, 120-175nm, Doug Kinnison                           !
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
        INTEGER, intent(in) :: nw
        INTEGER, intent(in) :: nz
        INTEGER, intent(inout) :: j
        REAL(rk), intent(in)    :: wl(:), wc(:)
        REAL(rk), intent(in)    :: tlev(:)
        REAL(rk), intent(in)    :: airden(:)
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg
        real(rk), optional, intent(out) :: sq(:,:)

!-----------------------------------------------------------------------------!
!     ... Local arrays                                                        !
!-----------------------------------------------------------------------------!
        
        integer, parameter :: kdata = 250
        integer n1, n2, n3, iskip, ierr
        integer i, iw, iz, n, idum, kk
        real(rk) :: x1    (kdata), x2   (kdata), x3(kdata)
        real(rk) :: y1    (kdata), y2   (kdata), y3(kdata)
        real(rk) :: qy(kz,kw)
        real(rk) :: so3  (kz,kw)
        real(rk) :: QY_O1D(nz,kw)
        real(rk), save :: yg298(kw), yg218(kw), xso3(kw)
        real(rk) :: tin(nz), yg(kw)
        real(rk) :: a1, a2, a3, w, t, kt, dum, q1, q2

        real(rk), parameter :: A(3) = (/0.8036_rk, 8.9061_rk, 0.1192_rk /)
        real(rk), parameter :: X(3) = (/ 304.225_rk, 314.957_rk, 310.737_rk/)
        real(rk), parameter :: om(3) = (/ 5.576_rk, 6.601_rk, 2.187_rk/)

        logical, save :: is_initialized = .false.
        integer :: chnl

        if (present(sq)) then
           sq = xnan
        end if

        if( initialize ) then
           if( .not. is_initialized ) then
              CALL readit
              is_initialized = .true.
           endif
        else
           chnl = xsqy_tab(j)%channel

           !-----------------------------------------------------------
           !     ... tin set to tlev
           !-----------------------------------------------------------
           tin(:nz) = tlev(:nz)
           !-------------------------------------------------------
           !     ... for hartley and huggins bands, use 
           !         temperature-dependent values from
           !         JPL06
           !-------------------------------------------------------
           !-------------------------------------------------------
           !     ... Cross Sections and Quantum Yields
           !-------------------------------------------------------
           do iw = 1, nw-1
              do iz = 1, nz

                 so3(iz,iw) = xso3(iw)

                 if ((wc(iw) .ge. 196.078_rk) .and. (wc(iw) .le. 342.5_rk)) then

                    if (tin(iz) .lt. 218._rk) then
                       so3(iz,iw) = yg218(iw)
                    endif
                    if ((tin(iz) .ge. 218._rk) .and. (tin(iz) .le. 298._rk)) then
                       so3(iz,iw) = yg218(iw)+(yg298(iw)-yg218(iw))/(298._rk-218._rk)* &
                            (tin(iz)-218._rk)
                    endif
                    if (tin(iz) .gt. 298._rk) then
                       so3(iz,iw) = yg298(iw)
                    endif
                 endif

              enddo
           enddo


           !------------------------------------------------------ 
           !     ... QY JPL06
           !         Valid from 306-328 nm
           !                    200-320 K
           !------------------------------------------------------ 
           do iz = 1, nz
              T = max(min(320.0_rk, tin(iz)),200._rk)

              do iw = 1, nw-1 

                 kt = 0.695_rk * T
                 q1 = 1._rk
                 q2 = exp(-825.518_rk/kT)

                 IF(wc(iw) .LE. 305._rk) THEN
                    QY_O1D (iz,iw) = 0.90_rk
                 ELSEIF((wc(iw) .GT. 305) .AND. (wc(iw) .LE. 328._rk)) THEN

                    QY_O1D(iz,iw)  = 0.0765_rk + &
                         a(1)*                (q1/(q1+q2))*EXP(-((x(1)-wc(iw))/om(1))**4)+ &
                         a(2)*(T/300._rk)**2 *(q2/(q1+q2))*EXP(-((x(2)-wc(iw))/om(2))**2)+ &
                         a(3)*(T/300._rk)**1.5_rk         *EXP(-((x(3)-wc(iw))/om(3))**2)

                 ELSEIF(wc(iw) .GT. 328._rk .AND. wc(iw) .LE. 345._rk) THEN
                    QY_O1D(iz,iw) = 0.08_rk
                 ELSEIF(wc(iw) .GT. 340._rk) THEN
                    QY_O1D(iz,iw) = 0._rk
                 ENDIF

                 QY_O1D(iz,iw) = min(qy_O1D(iz,iw),1.0_rk)
              enddo

           enddo
           !------------------------------------------------------
           !     ... derive the cross section*qy
           !------------------------------------------------------

           if (chnl == 1) then
              qy(:nz,:nw-1) = qy_O1D(:nz,:nw-1)
           else
              qy(:nz,:nw-1) = 1.0_rk - qy_O1D(:nz,:nw-1)
           end if

           do iz = 1, nz
              do iw = 1, nw-1
                 sq(iz,iw) = qy(iz,iw)*so3(iz,iw)
              enddo
           enddo
        end if

      contains

        subroutine readit
          !-----------------------------------------------------------
          !     Ref: Atkinson, ultraviolet solar radiation related to 
          !          mesopheric procesess, 149-159, in fiocco, g. (ed.), 
          !          mesospheric models and related exp., d. reidel, 
          !          dordrecht, 1971. 
          !
          !          120 nm through 200 nm
          !-----------------------------------------------------------
          x1 = xnan
          y1 = xnan
          n = xsqy_tab(j)%filespec%nread(1)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
               rd_cnt=n, &
               x=x1,y=y1 )
          call add_pnts_inter2(x1,y1,yg, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          do iw = 1, nw-1
             xso3(iw) = yg(iw)
          enddo

          !-------------------------------------------------------
          !     ... REF: JPL06 218K
          !         from 196.078 to 342.5 nm
          !-------------------------------------------------------
          x1 = xnan
          y1 = xnan
          y2 = xnan
          n = xsqy_tab(j)%filespec%nread(2)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(2), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(2), &
               rd_cnt=n, &
               x=x1,y=y1, y1=y2)
          x2(:n) = 0.5_rk*(x1(:n) + y1(:n))
          call add_pnts_inter2(x2,y2,yg218, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          !-------------------------------------------------------
          !     ... REF: JPL06 293-298K
          !         from 185.185 to 827.500 nm
          !-------------------------------------------------------
          x1 = xnan
          y1 = xnan
          y3 = xnan
          n = xsqy_tab(j)%filespec%nread(3)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(3), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(3), &
               rd_cnt=n, &
               x=x1,y=y1, y1=y3)
          x3(:n) = 0.5_rk*(x1(:n) + y1(:n))         
          call add_pnts_inter2(x3,y3,yg298, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          do iw = 1, nw-1
             if (wc(iw) .ge. 184.0_rk) then
                xso3(iw) = yg298(iw)
             endif
          enddo

       end subroutine readit

      end subroutine XSQY_O3

      subroutine XSQY_NO3(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide the product (absorptioon cross section) x (quantum yield) for     !
!   both channels of no3 photolysis:                                          !
!           (a) NO3 + hv -> NO2 + O(3P)                                       !
!           (b) NO3 + hv -> NO + O2                                           !
!   cross section and quantum yield consistent with JPL06                     !
!-----------------------------------------------------------------------------!
!   parameters:                                                               !
!   nw     - integer, number of specified intervals + 1 in working        (i) !
!            wavelength grid                                                  !
!   wl     - real, vector of lower limits of wavelength intervals in      (i) !
!            working wavelength grid                                          !
!   wc     - real, vector of center points of wavelength intervals in     (i) !
!            working wavelength grid                                          !
!   nz     - integer, number of altitude levels in working altitude grid  (i) !
!   tlev   - real, temperature (k) at each specified altitude level       (i) !
!   airlev - real, air density (molec/cc) at each altitude level          (i) !
!   j      - integer, counter for number of weighting functions defined  (io) !
!   sq     - real, cross section x quantum yield (cm^2) for each          (o) !
!            photolysis reaction defined, at each defined wavelength and      !
!            at each defined altitude level                                   !
!   jlabel - character*60, string identifier for each photolysis reaction (o) !
!            defined                                                          !
!-----------------------------------------------------------------------------!
!   edit history:                                                             !
!   02/04/08 Doug Kinnison                                                    !
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
        INTEGER, intent(in) :: nw
        INTEGER, intent(in) :: nz
        INTEGER, intent(inout) :: j
        REAL(rk), intent(in)    :: wl(:), wc(:)
        REAL(rk), intent(in)    :: tlev(:)
        REAL(rk), intent(in)    :: airden(:)
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg
        real(rk), optional, intent(out) :: sq(:,:)

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
        integer kdata
        parameter(kdata=350)
        integer i, iw, iz, n, n1, idum, iskip, ierr

        real(rk) :: x1(kdata)
        real(rk) :: y1(kdata), xqy( kdata)
        real(rk) :: qyNO_298 (kdata), qyNO_230 (kdata), qyNO_190 (kdata)
        real(rk) :: qyNO2_298(kdata), qyNO2_230(kdata), qyNO2_190(kdata)
        real(rk), save :: ygNO_298 (kw),    ygNO_230 (kw),    ygNO_190 (kw)
        real(rk), save :: ygNO2_298(kw),    ygNO2_230(kw),    ygNO2_190(kw)
        real(rk), save :: yg(kw)
        real(rk) :: qyNO2(nz,kw),     qyNO(nz,kz)
        real(rk) :: tin(nz)
        real(rk) :: qy
        logical, save :: is_initialized = .false.

        if (present(sq)) then
           sq = xnan
        end if

        if( initialize ) then
           if( .not. is_initialized ) then
              CALL readit
              is_initialized = .true.
           endif
        else

           !----------------------------------------------
           !     ... tin set to tlev
           !---------------------------------------------
           tin(:) = tlev(:)


           do iw = 1, nw-1

              do iz = 1, nz

                 IF (wc(iw) .GE. 585._rk) THEN
                    !... NO2 + O 
                    qyNO2(iz,iw) = ygNO2_190(iw)
                    if ((tin(iz) .GE. 190._rk) .AND. (tin(iz) .le. 230._rk)) then
                       qyNO2(iz,iw) = ygNO2_190(iw) + &
                            (ygNO2_230(iw)-ygNO2_190(iw))/(230._rk-190._rk) * &
                            (tin(iz)-190._rk)
                    endif
                    if ((tin(iz) .GT. 230._rk) .AND. (tin(iz) .le. 298._rk)) then
                       qyNO2(iz,iw) = ygNO2_230(iw) + &
                            (ygNO2_298(iw)-ygNO2_230(iw))/(298._rk-230._rk) * &
                            (tin(iz)-230._rk)
                    endif
                    if (tin(iz) .GT. 298._rk) then
                       qyNO2(iz,iw) = ygNO2_298(iw)
                    endif
                    !... NO + O2
                    qyNO(iz,iw) = ygNO_190(iw)
                    if ((tin(iz) .GE. 190._rk) .AND. (tin(iz) .le. 230._rk)) then
                       qyNO(iz,iw) = ygNO_190(iw) + &
                            (ygNO_230(iw)-ygNO_190(iw))/(230._rk-190._rk) * &
                            (tin(iz)-190._rk)
                    endif
                    if ((tin(iz) .GT. 230._rk) .AND. (tin(iz) .le. 298._rk)) then
                       qyNO(iz,iw) = ygNO_230(iw) + &
                            (ygNO_298(iw)-ygNO_230(iw))/(298._rk-230._rk) * &
                            (tin(iz)-230._rk)
                    endif
                    if (tin(iz) .GT. 298._rk) then
                       qyNO(iz,iw) = ygNO_298(iw)
                    endif

                 ELSE
                    qyNO(iz,iw) = 0._rk
                    qyNO2(iz,iw)= 1._rk
                 ENDIF
              enddo
           enddo
           if (xsqy_tab(j)%channel == 1) then
              ! 'NO3 + hv -> NO2 + O(3P)'
              do iw = 1, nw-1
                 do iz = 1, nz
                    sq(iz,iw) = qyNO(iz,iw) * yg(iw)

                 enddo
              enddo
           else
              ! 'NO3 + hv  -> NO + O2'
              do iw = 1, nw-1
                 do iz = 1, nz
                    sq(iz,iw) = qyNO2(iz,iw) * yg(iw)
                 enddo
              enddo
           end if
        end if

      contains

        subroutine readit
          x1 = xnan
          y1 = xnan
          n = xsqy_tab(j)%filespec%nread(1)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
               rd_cnt=n, &
               x=x1,y=y1 )

          call add_pnts_inter2(x1,y1,yg, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          xqy = xnan
          qyNO_298 = xnan
          qyNO_230 = xnan
          qyNO_190 = xnan
          qyNO2_298 = xnan
          qyNO2_230 = xnan
          qyNO2_190 = xnan
          n = xsqy_tab(j)%filespec%nread(2)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(2), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(2), &
               rd_cnt=n, &
               x=xqy, y=qyNO_298, y1=qyNO_230, y2=qyNO_190, &
               y3=qyNO2_298, y4=qyNO2_230, y5=qyNO2_190 )

          qyNO_298 = qyNO_298 * xsqy_tab(j)%filespec%xfac(2)
          qyNO_230 = qyNO_230 * xsqy_tab(j)%filespec%xfac(2)
          qyNO_190 = qyNO_190 * xsqy_tab(j)%filespec%xfac(2)
          qyNO2_298 = qyNO2_298 * xsqy_tab(j)%filespec%xfac(2)
          qyNO2_230 = qyNO2_230 * xsqy_tab(j)%filespec%xfac(2)
          qyNO2_190 = qyNO2_190 * xsqy_tab(j)%filespec%xfac(2)

          x1(:n) = xqy(:n)
          y1(:n) = qyNO2_298(:n)
          call add_pnts_inter2(x1,y1,ygNO2_298, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          y1(:n) = qyNO2_230(:n)
          call add_pnts_inter2(x1,y1,ygNO2_230, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          y1(:n) = qyNO2_190(:n)
          call add_pnts_inter2(x1,y1,ygNO2_190, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          y1(:n) = qyNO_298(:n)
          call add_pnts_inter2(x1,y1,ygNO_298, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          y1(:n) = qyNO_230(:n)
          call add_pnts_inter2(x1,y1,ygNO_230, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

          y1(:n) = qyNO_190(:n)
          call add_pnts_inter2(x1,y1,ygNO_190, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

        end subroutine readit

      end subroutine XSQY_NO3


      subroutine XSQY_N2O5(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield):                        !
!           N2O5 + hv -> NO3 + NO2                                            !
!           N2O5 + hv -. NO3 + NO + O                                         !
!   cross section: JPL06                                                      !
!-----------------------------------------------------------------------------!
!   parameters:                                                               !
!   nw     - integer, number of specified intervals + 1 in working        (i) !
!            wavelength grid                                                  !
!   wl     - real, vector of lower limits of wavelength intervals in      (i) !
!            working wavelength grid                                          !
!   wc     - real, vector of center points of wavelength intervals in     (i) !
!            working wavelength grid                                          !
!   nz     - integer, number of altitude levels in working altitude grid  (i) !
!   tlev   - real, temperature (k) at each specified altitude level       (i) !
!   airlev - real, air density (molec/cc) at each altitude level          (i) !
!   j      - integer, counter for number of weighting functions defined  (io) !
!   sq     - real, cross section x quantum yield (cm^2) for each          (o) !
!            photolysis reaction defined, at each defined wavelength and      !
!            at each defined altitude level                                   !
!   jlabel - character*60, string identifier for each photolysis reaction (o) !
!            defined                                                          !
!-----------------------------------------------------------------------------!
!   edit history:                                                             !
!   01/17/08  Doug Kinnison                                                   !
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
        INTEGER, intent(in) :: nw
        INTEGER, intent(in) :: nz
        INTEGER, intent(inout) :: j
        REAL(rk), intent(in)    :: wl(:), wc(:)
        REAL(rk), intent(in)    :: tlev(:)
        REAL(rk), intent(in)    :: airden(:)
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg
        real(rk), optional, intent(out) :: sq(:,:)


        !-----------------------------------------------------------------------------!
        !     ... local                                                               !
        !-----------------------------------------------------------------------------!
        integer, parameter :: kdata=300
        real(rk) :: x1   (kdata)
        real(rk) :: y1   (kdata)
        real(rk) :: wctmp(kdata)
        real(rk) :: wcb  (kdata)
        real(rk) :: ytmp (nz,kdata)
        real(rk) :: ycomb(nz,kdata)
        real(rk) :: ytd  (nz,kw)
        real(rk) :: yg   (kw)
        real(rk) :: yg1  (kw)
        real(rk) :: qy_O3p
        real(rk) :: tin(nz)
        real(rk) :: XS_harwood (nz,16)
        real(rk), save :: Xin(kdata) = -huge(1._rk)
        real(rk), save :: Yin(kdata) = -huge(1._rk)

        integer i, iw, n, idum, ierr
        integer iz, icnt, iwc, n1, chnl

        real(rk), parameter :: ww(16) = (/ &
             260.0_rk,   270.0_rk,   280.0_rk,    290.0_rk,   300.0_rk, &
             310.0_rk,   320.0_rk,   330.0_rk,    340.0_rk,   350.0_rk, &
             360.0_rk,   370.0_rk,   380.0_rk,    390.0_rk,   400.0_rk, &
             410.0_rk /)

        real(rk), parameter :: aa(16) = (/ &
             -18.27_rk,  -18.42_rk,  -18.59_rk,   -18.72_rk,  -18.84_rk, &
             -18.90_rk,  -18.93_rk,  -18.87_rk,   -18.77_rk,  -18.71_rk, &
             -18.31_rk,  -18.14_rk,  -18.01_rk,   -18.42_rk,  -18.59_rk, &
             -18.13_rk /)

        real(rk), parameter :: bb(16) = (/ &
             -0.091_rk,  -0.104_rk,  -0.112_rk,   -0.135_rk,  -0.170_rk,  &
             -0.226_rk,  -0.294_rk,  -0.388_rk,   -0.492_rk,  -0.583_rk,  &
             -0.770_rk,  -0.885_rk,  -0.992_rk,   -0.949_rk,  -0.966_rk, &
             -1.160_rk /)
        LOGICAL, save :: is_initialized = .false.

        if (present(sq)) then
           sq = xnan
        end if

        if( initialize ) then
           if( .not. is_initialized ) then
              CALL readit
              is_initialized = .true.
           endif
        else

           n = xsqy_tab(j)%filespec%nread(1)

           !----------------------------------------------------
           !     ... tin set to tlev
           !----------------------------------------------------
           tin(:) = tlev(:)

           !----------------------------------------------------
           !     ... Calculate the T-dep XS (233-295K)
           !         and 260-410nm
           !----------------------------------------------------
           do iw = 1, 16
              do iz = 1, nz
                 IF (tin(iz) .LT. 200.0_rk) THEN
                    XS_harwood(iz,iw) = 10**(aa(iw) + (1000._rk*bb(iw)/200.0_rk))
                 ENDIF
                 IF ((tin(iz) .GE. 200.0_rk) .AND. (tin(iz) .LE. 295._rk)) THEN
                    XS_harwood(iz,iw) = 10**(aa(iw) + (1000._rk*bb(iw)/tin(iz)))
                 ENDIF
                 IF (tin(iz) .GT. 295.0_rk) THEN
                    XS_harwood(iz,iw) = 10**(aa(iw) + (1000._rk*bb(iw)/295.0_rk))
                 ENDIF
              enddo
           enddo

           !     ... Combine cross sections
           do iz = 1, nz
              icnt = 1

              !     ... < 260 nm
              do i = 1, n
                 IF (xin(i) .LT. 260._rk) THEN
                    ycomb(iz,icnt) = yin(i)
                    wcb  (icnt)    = xin(i)
                    icnt = icnt + 1
                 ENDIF
              enddo
              !     ... 260-410 nm
              do i = 1, 16
                 ycomb(iz,icnt) = (xs_harwood(iz,i))
                 wcb  (icnt)    =  ww(i)
                 icnt = icnt+1
              enddo
              !     ... >410 nm
              do i = 1, n
                 IF (xin(i) .GT. 410._rk) THEN
                    ycomb(iz,icnt) = yin(i)
                    wcb  (icnt)    = xin(i)
                    icnt = icnt+1
                 ENDIF
              enddo
           enddo

           !     ... Interpolate to TUV grid 
           do iz = 1, nz
              n1 = icnt-1
              y1 = ycomb(iz,:)
              x1 = wcb
              call add_pnts_inter2(x1,y1,yg1, kdata, n1, &
                             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

              ytd(iz,:) = yg1(:)

           enddo

           !-------------------------------------------------------
           !     ... quantum yield (JPL06)
           !-------------------------------------------------------

           chnl = xsqy_tab(j)%channel
           do iw = 1, nw-1

              if (wc(iw) .GE. 300.0_rk) THEN 
                 qy_O3p = 0.0_rk
                 if (chnl.eq.2) then ! 'N2O5 + hv -> NO3 + NO2'
                    do iz = 1, nz
                       sq(iz,iw) = 1.0_rk * ytd(iz,iw)
                    enddo
                 else                ! 'N2O5 + hv -> NO3 + NO + O'
                    do iz = 1, nz
                       sq(iz,iw) = qy_O3p
                    enddo
                 endif
              endif

              if (wc(iw) .LT. 300.0_rk) THEN
                 qy_O3p = min( 1._rk, 3.832441_rk - 0.012809638_rk * wc(iw) )
                 qy_O3p = max( 0._rk, qy_O3p )
                 if (chnl.eq.2) then ! 'N2O5 + hv -> NO3 + NO2'
                    do iz = 1, nz
                       sq(iz,iw) = (1.0_rk-qy_O3p)*ytd(iz,iw)

                    enddo
                 else                ! 'N2O5 + hv -> NO3 + NO + O'
                    do iz = 1, nz
                       sq(iz,iw) = qy_O3p *ytd(iz,iw)
                    enddo
                 endif
              endif

           enddo

        endif

      contains
        subroutine readit
          n = xsqy_tab(j)%filespec%nread(1)
          CALL base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
               rd_cnt=n, &
               x=Xin,y=Yin )
        end subroutine readit
      end subroutine  XSQY_N2O5

      subroutine XSQY_CH3CHO(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!---------------------------------------------------------------------------!
!  PURPOSE:                                                                 !
!  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: !
!      (a)  CH3CHO + hv -> CH3 + HCO                                        !
!      (b)  CH3CHO + hv -> CH4 + CO                                         !
!      (c)  CH3CHO + hv -> CH3CO + H                                        !
!  Cross section:  Choice between                                           !
!                   (1) IUPAC 97 data, from Martinez et al.                 !
!                   (2) Calvert and Pitts                                   !
!                   (3) Martinez et al., Table 1 scanned from paper         !
!                   (4) KFA tabulations                                     !
!  Quantum yields: Choice between                                           !
!                   (1) IUPAC 97, pressure correction using Horowith and    !
!                                 Calvert, 1982                             !
!                   (2) NCAR data file, from Moortgat, 1986                 !
!---------------------------------------------------------------------------!
!  PARAMETERS:                                                              !
!  NW     - INTEGER, number of specified intervals + 1 in working        (I)!
!           wavelength grid                                                 !
!  WL     - REAL, vector of lower limits of wavelength intervals in      (I)!
!           working wavelength grid                                         !
!  WC     - REAL, vector of center points of wavelength intervals in     (I)!
!           working wavelength grid                                         !
!  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)!
!  TLEV   - REAL, temperature (K) at each specified altitude level       (I)!
!  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)!
!  J      - INTEGER, counter for number of weighting functions defined  (IO)!
!  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)!
!           photolysis reaction defined, at each defined wavelength and     !
!           at each defined altitude level                                  !
!  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)!
!           defined                                                         !
!---------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
      integer kdata
      parameter (kdata=150)

      integer  i, n, n1, n2, ierr, iz, iw, idum
      real(rk)  x1(kdata), x2(kdata)
      real(rk)  y1(kdata), y2(kdata)
      real(rk), save :: yg(kw) = -huge(1._rk)
      real(rk), save :: yg1(kw) = -huge(1._rk)
      real(rk), save :: yg2(kw) = -huge(1._rk)
      real(rk), save :: yg3(kw) = -huge(1._rk)
      real(rk), save :: yg4(kw) = -huge(1._rk)
      real(rk)  qy1, qy2, qy3
      real(rk)  sig

      if (present(sq)) then
         sq = xnan
      end if
      
!----------------------------------------------------
!... CH3CHO photolysis
!       1:  CH3 + HCO
!       2:  CH4 + CO
!       3:  CH3CO + H
!----------------------------------------------------
!      j = j+1
!      jlabel(j) = 'ch3cho -> ch3 + hco'	
!      j = j+1
!      jlabel(j) = 'ch3cho -> ch4 + co'	
!      j = j+1
!      jlabel(j) = 'ch3cho -> ch3co + h'
!----------------------------------------------------

      if( initialize ) then
         CALL readit
      else
         !----------------------------------------------------
         !...  combine XS*QY 
         !----------------------------------------------------
         DO iw = 1, nw - 1
            DO i = 1, nz

               sig = yg(iw)
               qy1 = yg1(iw)
               qy2 = yg2(iw)
               qy3 = yg3(iw)

               !... Pressure correction for channel 1, CH3 + CHO
               !     based on Horowitz and Calvert 1982.

               qy1 = qy1 * (1._rk + yg4(iw))/(1._rk + yg4(iw)*airden(i)/2.465E19_rk)
               qy1 = MIN(1._rk, qy1)
               qy1 = MAX(0._rk, qy1)

               sq(i,iw) = (sig*qy1)+(sig*qy2) + (sig * qy3)

            ENDDO
         ENDDO
      endif

    contains
      
      subroutine readit
        x1 = xnan
        y1 = xnan
        n = xsqy_tab(j)%filespec%nread(1)
        CALL base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
             errmsg=errmsg, errflg=errflg, &
             skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
             rd_cnt=n, &
             x=x1,y=y1 )
        call add_pnts_inter2(x1,y1,yg, kdata, n, &
             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

        x1 = xnan
        y1 = xnan
        y2 = xnan

        n = xsqy_tab(j)%filespec%nread(2)
        CALL base_read( filespec=xsqy_tab(j)%filespec%filename(2), &
             errmsg=errmsg, errflg=errflg, &
             skip_cnt=xsqy_tab(j)%filespec%nskip(2), &
             rd_cnt=n, &
             x=x1,y=y2, y1=y1 )
        call add_pnts_inter2(x1,y1,yg1, kdata, n, &
             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

        call add_pnts_inter2(x1,y2,yg2, kdata, n, &
             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
        do iw = 1, nw-1
           yg3(iw) = 0._rk
        enddo

        x1 = xnan
        y1 = xnan
        n = xsqy_tab(j)%filespec%nread(3)
        CALL base_read( filespec=xsqy_tab(j)%filespec%filename(3), &
             errmsg=errmsg, errflg=errflg, &
             skip_cnt=xsqy_tab(j)%filespec%nskip(3), &
             rd_cnt=n, &
             x=x1,y=x2, y1=y2, y2=y1 )        
        call add_pnts_inter2(x1,y1,yg4, kdata, n, &
             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

      end subroutine readit


    end subroutine XSQY_CH3CHO


          SUBROUTINE XSQY_MGLY(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!---------------------------------------------------------------------------!
!  PURPOSE:                                                                 !
!  Provide the product (cross section) x (quantum yield) for CH3COCHO       !
!  photolysis:                                                              !
!           MGLY (CH3COCHO) + hv -> CH3CO + HCO                             !
!                                                                           !
!  Cross section:                                                           !
!                Average at 1 nm of Staffelbach et al., 1995, and           !
!                      Meller et al., 1991                                  !
!  Quantum yield:                                                           !
!                 Chen, Y., W. Wang, and L. Zhu, Wavelength-dependent       !
!                 photolysis of methylglyoxal in the 290-440 nm region,     !
!                 J Phys Chem A, 104, 11126-11131, 2000.                    !
!---------------------------------------------------------------------------!
!  PARAMETERS:                                                              !
!  NW     - INTEGER, number of specified intervals + 1 in working        (I)!
!           wavelength grid                                                 !
!  WL     - REAL, vector of lower limits of wavelength intervals in      (I)!
!           working wavelength grid                                         !
!  WC     - REAL, vector of center points of wavelength intervals in     (I)!
!           working wavelength grid                                         !
!  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)!
!  TLEV   - REAL, temperature (K) at each specified altitude level       (I)!
!  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)!
!  J      - INTEGER, counter for number of weighting functions defined  (IO)!
!  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)!
!           photolysis reaction defined, at each defined wavelength and     !
!           at each defined altitude level                                  !
!  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)!
!           defined                                                         !
!---------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

!---------------------------------------------------------------------------!
!     ... local                                                             !
!---------------------------------------------------------------------------!
      integer kdata
      parameter (kdata=500)

      integer  i, n, n1, n2, ierr, iw, iz
      real(rk) :: x1(kdata), y1(kdata)
      real(rk),save :: yg(kw)=-huge(1._rk)
      real(rk) :: qy
      real(rk) :: sig
      real(rk) :: phi0, kq

      if (present(sq)) then
         sq = xnan
      end if
      if( initialize ) then
         CALL readit
      else

         !-------------------------------------------------------
         !     ... combine qy * xs
         !-------------------------------------------------------
         DO iw = 1, nw - 1

            sig = yg(iw)

            DO i = 1, nz

               !----------------------------------------
               !             quantum yields:
               !             zero pressure yield:
               !             1.0 for wc < 380 nm
               !             0.0 for wc > 440 nm
               !             linear in between:
               !----------------------------------------
               phi0 = 1._rk - (wc(iw) - 380._rk)/60._rk
               phi0 = MIN(phi0,1._rk)
               phi0 = MAX(phi0,0._rk)

               !----------------------------------------------------------
               !              Pressure correction: 
               !              quenching coefficient, torr-1
               !              in air, Koch and Moortgat:
               !----------------------------------------------------------
               kq = 1.36e8_rk * EXP(-8793_rk/wc(iw))

               !----------------------------------------------------------
               !              In N2, Chen et al:
               !----------------------------------------------------------
               !              kq = 1.93e4 * EXP(-5639/wc(iw))

               IF(phi0 .GT. 0._rk) THEN
                  IF (wc(iw) .GE. 380._rk .AND. wc(iw) .LE. 440._rk) THEN
                     qy = phi0 / (phi0 + kq * airden(i) * 760._rk/2.456E19_rk)
                  ELSE
                     qy = phi0
                  ENDIF
               ELSE
                  qy = 0._rk
               ENDIF

               sq(i,iw) = sig * qy

            ENDDO
         ENDDO


      endif


      contains
      
        subroutine readit
          x1 = xnan
          y1 = xnan
          n = xsqy_tab(j)%filespec%nread(1)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
               rd_cnt=n, &
               x=x1,y=y1 )
          call add_pnts_inter2(x1,y1,yg, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
        end subroutine readit

      end subroutine XSQY_MGLY

      SUBROUTINE XSQY_ACETONE(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!---------------------------------------------------------------------------!
!  PURPOSE:                                                                 !
!  Provide product (cross section) x (quantum yield) for CH3COCH3 photolysis!
!          CH3COCH3 + hv -> Products                                        !
!                                                                           !
!  Cross section:  Choice between                                           !
!                   (1) Calvert and Pitts                                   !
!                   (2) Martinez et al., 1991, alson in IUPAC 97            !
!                   (3) NOAA, 1998, unpublished as of 01/98                 !
!  Quantum yield:  Choice between                                           !
!                   (1) Gardiner et al, 1984                                !
!                   (2) IUPAC 97                                            !
!                   (3) McKeen et al., 1997                                 !
!---------------------------------------------------------------------------!
!  PARAMETERS:                                                              !
!  NW     - INTEGER, number of specified intervals + 1 in working        (I)!
!           wavelength grid                                                 !
!  WL     - REAL, vector of lower limits of wavelength intervals in      (I)!
!           working wavelength grid                                         !
!  WC     - REAL, vector of center points of wavelength intervals in     (I)!
!           working wavelength grid                                         !
!  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)!
!  TLEV   - REAL, temperature (K) at each specified altitude level       (I)!
!  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)!
!  J      - INTEGER, counter for number of weighting functions defined  (IO)!
!  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)!
!           photolysis reaction defined, at each defined wavelength and     !
!           at each defined altitude level                                  !
!  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)!
!           defined                                                         !
!---------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
      integer  kdata
      parameter (kdata=150)

      integer  i, n, n1, n2, n3, iw, ierr, iz, idum
      real(rk),save :: x1(kdata)=-huge(1._rk)
      real(rk),save :: y1(kdata)=-huge(1._rk)
      real(rk),save :: A(kdata)=-huge(1._rk)
      real(rk),save :: B(kdata)=-huge(1._rk)
      real(rk),save :: C(kdata)=-huge(1._rk)
      real(rk) ::  x2(kdata), y2(kdata)
      real(rk) ::  xs(nz,kdata), sig(nz,kw)
      real(rk) ::  yg(kw), yg1(kw), yg2(kw), yg3(kw)
      real(rk) ::  tin(nz), AD(nz)
      real(rk) ::  qytot(kw), qyCO(kw), qyCH3CO(kw)
      real(rk) ::  AA0, a0, b0
      real(rk) ::  AA1, a1, b1, t, qy
      real(rk) ::  AA2, AA3, AA4, a2, b2, a3, b3, c3, a4, b4

      if (present(sq)) then
         sq = xnan
      end if
      
!---------------------------------------------
!     ... CH3COCH3 photodissociation
!---------------------------------------------

      if( initialize ) then
         call readit
      else
!---------------------------------------------
!     ... tin set to tlev
!---------------------------------------------
         tin(:) = tlev(:)
         AD (:) = airden(:)

         n1 = xsqy_tab(j)%filespec%nread(1)
!---------------------------------------------
!     ... Derive XS at given temperature
!---------------------------------------------
    
         do iz = 1, nz

            do iw = 1, n1

               if ((tin(iz) .GE. 235._rk) .AND. (tin(iz) .LE. 298._rk)) Then
                  xs(iz,iw) = y1(iw) *( 1 + (A(iw)*tin(iz)) + &
                       (B(iw)*tin(iz)**2)  + &
                       (C(iw)*tin(iz)**3) )

               endif

               if (tin(iz) .LT. 235._rk) then
                  xs(iz,iw) = y1(iw) *( 1 + (A(iw)*235._rk) + &
                       (B(iw)*(235._rk)**2)  + &
                       (C(iw)*(235._rk)**3) )

               endif

               if (tin(iz) .GT. 298._rk) then
                  xs(iz,iw) = y1(iw) *( 1 + (A(iw)*298._rk) + &
                       (B(iw)*(298._rk)**2)  + &
                       (C(iw)*(298._rk)**3) )

               endif

            enddo

            n     = n1

            x2(:) = x1(:)
            y2(:) = xs(iz,:)

            call add_pnts_inter2(x2,y2,yg, kdata, n, &
                 nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)


            sig(iz,:) = yg(:)

         enddo

!---------------------------------------------
!     ... quantum yield JPL06
!---------------------------------------------
         DO iz = 1, nz

            T = min(tin(iz), 295._rk)
            T = max(T, 218._rk)

            DO iw = 1, nw-1

               qyCO(iw) = 0._rk
               qyCH3CO(iw) =  0._rk

               IF ((wc(iw) .GE. 279._rk).AND.(wc(iw) .LT. 327._rk) ) THEN

                  a0 = 0.350_rk* (T/295._rk)**(-1.28_rk)
                  b0 = 0.068_rk* (T/295._rk)**(-2.65_rk)
                  AA0 = (a0 / (1._rk-a0))* exp(b0*(wc(iw)-248._rk))
                  qyCO(iw) = 1._rk / (1._rk + AA0)
               ENDIF

               IF ((wc(iw) .GE. 279._rk).AND.(wc(iw) .LT. 302._rk)) THEN

                  a1 = 1.6e-19_rk* (T/295._rk)**(-2.38_rk) 
                  b1 = 0.55e-3_rk* (T/295._rk)**(-3.19_rk)
                  AA1 = a1* exp(-b1*((1e7_rk/wc(iw)) - 33113._rk))
                  qyCH3CO(iw) = (1._rk-qyCO(iw)) / (1._rk + AA1*AD(iz))

               ELSEIF ((wc(iw) .GE. 302._rk).AND.(wc(iw) .LE. 327.5_rk)) THEN

                  a2= 1.62e-17_rk* (T/295._rk)**(-10.03_rk)
                  b2= 1.79e-3_rk * (T/295._rk)**(-1.364_rk)
                  AA2= a2* exp(-b2*((1e7_rk/wc(iw))-30488._rk))

                  a3= 26.29_rk*   (T/295._rk)**(-6.59_rk)
                  b3= 5.72e-7_rk* (T/295._rk)**(-2.93_rk)
                  c3= 30006._rk*  (T/295._rk)**(-0.064_rk)
                  AA3= a3* exp(-b3*((1e7_rk/wc(iw))-c3)**2)

                  a4= 1.67e-15_rk* (T/295._rk)**(-7.25_rk)
                  b4= 2.08e-3_rk*  (T/295._rk)**(-1.16_rk)
                  AA4= a4* exp(-b4*((1e7_rk/wc(iw)) - 30488._rk))

                  qyCH3CO(iw) = ((1._rk + AA4*AD(iz) + AA3) / &
                       ((1._rk + AA2*AD(iz) + AA3)* &
                       (1._rk + AA4*AD(iz))))*(1-qyCO(iw))

               ENDIF

               qytot(iw) = qyCO(iw) + qyCH3CO(iw)

               if (wc(iw) .LT. 279._rk) then
                  qytot(iw) = 1.0_rk
               endif

               sq(iz,iw) = sig(iz,iw)*qytot(iw)

            ENDDO
         ENDDO
      endif      

      contains
      
        subroutine readit
          x1 = xnan
          y1 = xnan
          n = xsqy_tab(j)%filespec%nread(1)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
               rd_cnt=n, &
               x=x1,y=y1 )

          n = xsqy_tab(j)%filespec%nread(2)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(2), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(2), &
               rd_cnt=n, &
               x=x1,y=A,y1=B,y2=C )
          
          A(:n) = A(:n)*1.e-3_rk
          B(:n) = B(:n)*1.e-5_rk
          C(:n) = C(:n)*1.e-8_rk
        end subroutine readit

      end subroutine XSQY_ACETONE

      SUBROUTINE XSQY_PAN(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!---------------------------------------------------------------------------!
!  PURPOSE:                                                                 !
!  Provide product (cross section) x (quantum yield) for PAN photolysis:    !
!       PAN + hv -> Products                                                !
!                                                                           !
!  Cross section: from Talukdar et al., 1995                                !
!  Quantum yield: Assumed to be unity                                       !
!---------------------------------------------------------------------------!
!  PARAMETERS:                                                              !
!  NW     - INTEGER, number of specified intervals + 1 in working        (I)!
!           wavelength grid                                                 !
!  WL     - REAL, vector of lower limits of wavelength intervals in      (I)!
!           working wavelength grid                                         !
!  WC     - REAL, vector of center points of wavelength intervals in     (I)!
!           working wavelength grid                                         !
!  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)!
!  TLEV   - REAL, temperature (K) at each specified altitude level       (I)!
!  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)!
!  J      - INTEGER, counter for number of weighting functions defined  (IO)!
!  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)!
!           photolysis reaction defined, at each defined wavelength and     !
!           at each defined altitude level                                  !
!  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)!
!           defined                                                         !
!---------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

!---------------------------------------------------------------------------!
!     ... local                                                             !
!---------------------------------------------------------------------------!
      integer   kdata
      parameter (kdata=100)
      integer i, iw, n, n2, ierr, iz, idum
      real(rk) ::  x1(kdata), x2(kdata)
      real(rk) ::  y1(kdata), y2(kdata)
      real(rk),save ::  yg(kw)= -huge(1._rk)
      real(rk),save :: yg2(kw) = -huge(1._rk)
      real(rk) ::  tin(nz)
      real(rk) ::  qy, sig

      if (present(sq)) then
         sq = xnan
      end if

      if( initialize ) then
         call readit
      else
         !----------------------------------------------
         !     ... tin set to tlev
         !----------------------------------------------
         tin(:) = tlev(:)

         !----------------------------------------------
         !    ... Quantum yield
         !----------------------------------------------
         qy = 1.0_rk

         DO iw = 1, nw-1
            DO iz = 1, nz

               sig = yg(iw) * EXP(yg2(iw)*(tin(iz)-298._rk))

               sq(iz,iw) = qy * sig

            ENDDO
         ENDDO

      endif

      contains

        subroutine readit
          x1 = xnan
          y1 = xnan
          y2 = xnan
          n = xsqy_tab(j)%filespec%nread(1)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
               rd_cnt=n, &
               x=x1, y=y1, y1=y2 )
          y1(:n) = y1(:n) * 1.E-20_rk
          y2(:n) = y2(:n) * 1E-3_rk
          x2(:n) = x1(:n)

          call add_pnts_inter2(x1,y1,yg, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
          call add_pnts_inter2(x2,y2,yg2, kdata, n, &
               nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
        end subroutine readit

      end subroutine XSQY_PAN
            
      SUBROUTINE XSQY_H2O(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield) for hcl photolysis:     !
!           H2O + hv -> products                                              !
!   cross section: taken from three sources                                   !
!     1) JPL06 (jpl97-4), 175.5 - 189.3                                       !
!     2) Cantrell et al., grl, 24, 17, 2195-2198, 1997,  183.0 - 193.0 nm     !
!     3) Yoshino et al.,  chemical physics, 211 (1996) 387-391, 120.38-188.03 !
!                                                                             !
!   quantum yield: is unity between 175.5 and 189.3                           !
!-----------------------------------------------------------------------------!
!   parameters:                                                               !
!   nw     - integer, number of specified intervals + 1 in working        (i) !
!            wavelength grid                                                  !
!   wl     - real, vector of lower limits of wavelength intervals in      (i) !
!            working wavelength grid                                          !
!   wc     - real, vector of center points of wavelength intervals in     (i) !
!            working wavelength grid                                          !
!   nz     - integer, number of altitude levels in working altitude grid  (i) !
!   tlev   - real, temperature (k) at each specified altitude level       (i) !
!   airlev - real, air density (molec/cc) at each altitude level          (i) !
!   j      - integer, counter for number of weighting functions defined  (io) !
!   sq     - real, cross section x quantum yield (cm^2) for each          (o) !
!            photolysis reaction defined, at each defined wavelength and      !
!            at each defined altitude level                                   !
!   jlabel - character*60, string identifier for each photolysis reaction (o) !
!            defined                                                          !
!-----------------------------------------------------------------------------!
!   edit history:                                                             !
!   06/11/01   original, dek addition                                         !
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
      integer kdata
      parameter(kdata=7000)
      integer n, i, iw, n1, n2, n3, idum, ierr, iz
      real(rk) :: x1(kdata), y1(kdata)
      real(rk) :: x2(kdata), y2(kdata)
      real(rk) :: x3(kdata), y3(kdata)
      real(rk), save :: sqx(kw,3)=-huge(1._rk)
      real(rk) :: yg(kw), yg1(kw), yg2(kw), yg3(kw)
      real(rk) :: qy

      integer :: chnl

      LOGICAL, save :: is_initialized = .false.

      if (present(sq)) then
         sq = xnan
      end if

      if( initialize ) then
         if( .not. is_initialized ) then
            CALL readit
            is_initialized = .true.
         endif
      else

         chnl = xsqy_tab(j)%channel
         sq(:nw-1,1) = sqx(:nw-1,chnl)

      end if

    contains
      subroutine readit

        x1 = xnan
        y1 = xnan
        n = xsqy_tab(j)%filespec%nread(1)
        call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
             errmsg=errmsg, errflg=errflg, &
             skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
             rd_cnt=n, &
             x=x1,y=y1 )

        y1(1:n) = y1(1:n) * xsqy_tab(j)%filespec%xfac(1)

        call add_pnts_inter2(x1,y1,yg1, kdata, n, &
             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

        x2 = xnan
        y2 = xnan
        n = xsqy_tab(j)%filespec%nread(2)
        call base_read( filespec=xsqy_tab(j)%filespec%filename(2), &
             errmsg=errmsg, errflg=errflg, &
             skip_cnt=xsqy_tab(j)%filespec%nskip(2), &
             rd_cnt=n, &
             x=x2,y=y2 )

        y2(1:n) = y2(1:n) * xsqy_tab(j)%filespec%xfac(2)

        call add_pnts_inter2(x2,y2,yg2, kdata, n, &
             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

        x3 = xnan
        y3 = xnan
        n = xsqy_tab(j)%filespec%nread(3)
        call base_read( filespec=xsqy_tab(j)%filespec%filename(3), &
             errmsg=errmsg, errflg=errflg, &
             skip_cnt=xsqy_tab(j)%filespec%nskip(3), &
             rd_cnt=n, &
             x=x3,y=y3 )

        call add_pnts_inter2(x3,y3,yg3, kdata, n, &
             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)


        !--------------------------------------------------------
        !     ...combine data sets (i.e., Yoshino et al., 1996 
        !        and Cantrell et al., 1997)
        !--------------------------------------------------------    
        do i = 1, nw-1
           if (wc(i) .lt. 183.0_rk) then
              yg(i) = yg3(i)
           elseif (wc(i) .le. 194.0_rk) then
              yg(i) = yg2(i)
           else
              yg(i) = 0._rk
           endif
        enddo

        !------------------------------------------------------
        ! quantum yield assumed to be unity (jpl97-4)
        !------------------------------------------------------
        ! 105 to 145 nm
        ! (JPL 1997 which references Stief, L.J., W.A. 
        ! Payne, and R. B. Klemm, A flash
        ! photolysis-resonance fluoresence study of the 
        ! formation of O(1D) in the photolysis of water 
        ! and the reaction of O(1D) with H2, Ar, and He, 
        ! J. Chem. Phys., 62, 4000, 1975.)
        sqx = 0._rk
        do iw = 1, nw-1

           if (wc(iw) .le. 145.0_rk) then

              sqx(iw,1) = yg(iw) * 0.890_rk
              sqx(iw,2) = yg(iw) * 0.110_rk
              sqx(iw,3) = yg(iw) * 0.0_rk

           end if

           !     ... > 145nm
           !         JPL97
           if (wc(iw) .gt. 145.0_rk) then

              sqx(iw,1) = yg(iw) * 1.0_rk
              sqx(iw,2) = yg(iw) * 0.0_rk
              sqx(iw,3) = yg(iw) * 0.0_rk

           end if

        end do  ! end wavelength loop

        ! Overwrite Lyamn Alpha
        ! Slanger, T.G., and G. Black, Photodissociative 
        ! channels at 1216A for H2O, NH3 and CH4,
        ! J. Chem. Phys., 77, 2432, 1982.)

        sqx(la_ndx,1) = yg(la_ndx) * 0.780_rk
        sqx(la_ndx,2) = yg(la_ndx) * 0.100_rk
        sqx(la_ndx,3) = yg(la_ndx) * 0.120_rk

      end subroutine readit


      end subroutine XSQY_H2O     
      
      SUBROUTINE XSQY_HO2NO2(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product of (cross section) x (quantum yield) for hno4 photolysis  !
!    chnl 1)   HO2NO2 + hv -> OH +  NO3                                           !
!    chnl 2)   HO2NO2 + hv -> HO2 + NO2                                           !
!   cross sections and QY from JPL06                                          !
!-----------------------------------------------------------------------------!
!   parameters:                                                               !
!   nw     - integer, number of specified intervals + 1 in working        (i) !
!            wavelength grid                                                  !
!   wl     - real, vector of lower limits of wavelength intervals in      (i) !
!            working wavelength grid                                          !
!   wc     - real, vector of center points of wavelength intervals in     (i) !
!            working wavelength grid                                          !
!   nz     - integer, number of altitude levels in working altitude grid  (i) !
!   tlev   - real, temperature (k) at each specified altitude level       (i) !
!   airlev - real, air density (molec/cc) at each altitude level          (i) !
!   j      - integer, counter for number of weighting functions defined  (io) !
!   sq     - real, cross section x quantum yield (cm^2) for each          (o) !
!            photolysis reaction defined, at each defined wavelength and      !
!            at each defined altitude level                                   !
!   jlabel - character*60, string identifier for each photolysis reaction (o) !
!            defined                                                          !
!-----------------------------------------------------------------------------!
!   edit history:                                                             !
!   05/98  original, adapted from former jspec1 subroutine                    !
!   06/01  modified by doug kinnison                                          !
!   01/08  modified by Doug Kinnison                                          !
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
      integer kdata
      parameter(kdata=100)
      integer i, iw, iz, n, n1, idum, ierr, icnt
      real(rk),save :: x1   (kdata)=-huge(1._rk), x2(kdata)=-huge(1._rk)
      real(rk) :: wcb(kdata), x(kdata), y(kdata)
      real(rk),save :: y1   (kdata)=-huge(1._rk), aa(kdata)=-huge(1._rk), bb (kdata)=-huge(1._rk)
      real(rk) :: ytmp (nz,kdata), ycomb(nz,kdata)
      real(rk) :: ytd  (nz,kw), yg(kw)
      real(rk) :: Q(nz), tin(nz), t

      LOGICAL, save :: is_initialized = .false.
      integer :: chnl
      
      if (present(sq)) then
         sq = xnan
      end if

      if( initialize ) then
         if( .not. is_initialized ) then
            CALL readit
            is_initialized = .true.
         endif
      else

         chnl = xsqy_tab(j)%channel

         !----------------------------------------------
         !     ... tin set to tlev
         !----------------------------------------------
         tin(:) = tlev(:) 


         n = xsqy_tab(j)%filespec%nread(1)
         n1 = xsqy_tab(j)%filespec%nread(2)

         !----------------------------------------------
         !     ...Derive T-dep Burkholder et al., 2002.)
         !----------------------------------------------
         do iz = 1, nz
            do iw = 1, n1
               t           = MAX(280._rk,MIN(tin(iz),350._rk))
               Q(iz)       = 1 + exp(-988._rk/(0.69_rk*t))
               ytmp(iz,iw) = ( aa(iw)/Q(iz) + bb(iw)*(1-1/Q(iz)))*1e-20_rk
            enddo
         enddo
         !     ... Combine cross sections
         do iz = 1, nz
            icnt = 1

            !     ... < 280 nm
            !     ... x1(iw) goes from 190-350nm
            do iw = 1, n
               IF (x1(iw) .LT. 280._rk) THEN
                  ycomb(iz,icnt) = y1(iw)
                  wcb  (icnt)    = x1(iw)
                  icnt = icnt + 1
               ENDIF
            enddo
            !     ... 280-350 nm
            do iw = 1, n1
               ycomb(iz,icnt) = ytmp(iz,iw)
               wcb  (icnt)    = x2  (iw)
               icnt = icnt+1
            enddo
         enddo

         !     ... Interpolate Combine cross sections
         do iz = 1, nz
            n  = icnt-1
            y = ycomb(iz,:)
            x = wcb

            call add_pnts_inter2(x,y,yg, kdata, n, &
                 nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

            ytd(iz,:) = yg(:)

         enddo

         do iw = 1, nw - 1
            IF (wc(iw) .LT. 200.0_rk) THEN
               if (chnl==1) then
                  sq(:nz,iw) = 0.30_rk * ytd(:nz,iw)
               else if (chnl==2) then
                  sq(:nz,iw) = 0.70_rk * ytd(:nz,iw)
               end if
            ELSEIF (wc(iw) .GE. 200.0_rk) THEN
               if (chnl==1) then
                  sq(:nz,iw) = 0.20_rk * ytd(:nz,iw)
               else if (chnl==2) then
                  sq(:nz,iw) = 0.80_rk * ytd(:nz,iw)
               end if
            ENDIF

         enddo
      endif

      contains

        subroutine readit
          x1 = xnan
          y1 = xnan

          n = xsqy_tab(j)%filespec%nread(1)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
               rd_cnt=n, &
               x=x1, y=y1 )
          
          x2 = xnan
          aa = xnan
          bb = xnan
          n = xsqy_tab(j)%filespec%nread(2)
          call base_read( filespec=xsqy_tab(j)%filespec%filename(2), &
               errmsg=errmsg, errflg=errflg, &
               skip_cnt=xsqy_tab(j)%filespec%nskip(2), &
               rd_cnt=n, &
               x=x2, y=aa, y1=bb )

        end subroutine readit
        
      end subroutine XSQY_HO2NO2

      subroutine XSQY_NOp(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield):                        !
!           NO + hv = NOp + e                                                 !
!   cross section: JPL06                                                      !
!   quantum yield: is unity.                                                  !
!-----------------------------------------------------------------------------!
!   parameters:                                                               !
!   nw     - integer, number of specified intervals + 1 in working        (i) !
!            wavelength grid                                                  !
!   wl     - real, vector of lower limits of wavelength intervals in      (i) !
!            working wavelength grid                                          !
!   wc     - real, vector of center points of wavelength intervals in     (i) !
!            working wavelength grid                                          !
!   nz     - integer, number of altitude levels in working altitude grid  (i) !
!   tlev   - real, temperature (k) at each specified altitude level       (i) !
!   airlev - real, air density (molec/cc) at each altitude level          (i) !
!   j      - integer, counter for number of weighting functions defined  (io) !
!   sq     - real, cross section x quantum yield (cm^2) for each          (o) !
!            photolysis reaction defined, at each defined wavelength and      !
!            at each defined altitude level                                   !
!   jlabel - character*60, string identifier for each photolysis reaction (o) !
!            defined                                                          !
!-----------------------------------------------------------------------------!
!   edit history:                                                             !
!   01/16/08  Doug Kinnison                                                   !
!-----------------------------------------------------------------------------!


!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
        INTEGER, intent(in) :: nw
        INTEGER, intent(in) :: nz
        INTEGER, intent(inout) :: j
        REAL(rk), intent(in)    :: wl(:), wc(:)
        REAL(rk), intent(in)    :: tlev(:)
        REAL(rk), intent(in)    :: airden(:)
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg
        real(rk), optional, intent(out) :: sq(:,:)

        !----------------------------------------------------
        !  M. Nicolet, "Aeronomical Aspects of Mesopheric Photodissociation: Prosesses Resulting
        !  from the Solar H Lyman-Alpha Line", Planet. Space Sci. Vol. 33, No. 1, pp. 69-80, 1985.
        !  see eq. (27)
        !----------------------------------------------------
        if( .not. initialize ) then

           sq(:,:) = 0._rk
           sq(la_ndx,1) = 2.02e-18_rk ! lyman alpha band

        endif

      end subroutine XSQY_NOp

      subroutine XSQY_CH2BR2(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield) for ch2br2 photolysis:  !
!          CH2Br2 + hv -> 2Br                                                 !
!   cross section: from JPL06 recommendation                                  !
!   quantum yield: assumed to be unity                                        !
!-----------------------------------------------------------------------------!
!   parameters:                                                               !
!   nw     - integer, number of specified intervals + 1 in working        (i) !
!            wavelength grid                                                  !
!   wl     - real, vector of lower limits of wavelength intervals in      (i) !
!            working wavelength grid                                          !
!   wc     - real, vector of center points of wavelength intervals in     (i) !
!            working wavelength grid                                          !
!   nz     - integer, number of altitude levels in working altitude grid  (i) !
!   tlev   - real, temperature (k) at each specified altitude level       (i) !
!   airlev - real, air density (molec/cc) at each altitude level          (i) !
!   j      - integer, counter for number of weighting functions defined  (io) !
!   sq     - real, cross section x quantum yield (cm^2) for each          (o) !
!            photolysis reaction defined, at each defined wavelength and      !
!            at each defined altitude level                                   !
!   jlabel - character*60, string identifier for each photolysis reaction (o) !
!            defined                                                          !
!-----------------------------------------------------------------------------!
!   edit history:                                                             !
!   07/30/07  Doug Kinnison                                                   !
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      !-----------------------------------------------------------------------------!
      !     ... local                                                               !
      !-----------------------------------------------------------------------------!
      integer kdata
      parameter(kdata=300)
      integer i, iw, n, idum, nloop, n1
      integer ierr, iz, iwc, icnt
      real(rk) :: x1   (kdata),   y1   (kdata)
      real(rk), save :: xin  (kdata),   yin  (kdata)
      real(rk) :: wctmp(kdata),   wcb  (kdata)
      real(rk) :: ytmp (nz,kdata),ycomb(nz,kdata)
      real(rk) :: yg1  (kw),      ytd  (nz,kw)
      real(rk) :: qy
      real(rk) :: tin(nz)

      real(rk), parameter :: AA(5) = (/ &
           -70.211776_rk , &
           1.940326e-1_rk , &
           2.726152e-3_rk , &
           -1.695472e-5_rk , &
           2.500066e-8_rk /)

      real(rk), parameter ::  BB(5) = (/ &
           2.899280_rk , &
           -4.327724e-2_rk , &
           2.391599e-4_rk , &
           -5.807506e-7_rk , &
           5.244883e-10_rk /)

      real(rk), parameter ::  lp(5) = (/ &
           0.0_rk , &
           1.0_rk , &
           2.0_rk , &
           3.0_rk , &
           4.0_rk /)

      n = xsqy_tab(j)%filespec%nread(1)

      if (present(sq)) then
         sq = xnan
      end if
      if( initialize ) then
         CALL readit
      else

         !----------------------------------------------
         !     ... set tin to tlev
         !----------------------------------------------
         tin(:)   = tlev(:)

         !----------------------------------------------
         !    Derive temperature dependence 
         !----------------------------------------------
         !    Temperature dependence good between 
         !      210-300K and 210 nm-290 nm
         !----------------------------------------------
         iwc      = 1
         ytmp(:,:)= 0.0_rk

         do iw = 1, nw-1

            IF ((wc(iw) .GE. 210._rk) .AND. (wc(iw) .LE.290._rk)) THEN

               do iz = 1, nz

                  IF (tin(iz) .LT. 210._rk) THEN
                     do nloop = 1, 5
                        ytmp(iz,iwc) = ytmp(iz,iwc) &
                             +  AA(nloop)* (wc(iw)**lp(nloop)) &
                             + (210.0_rk-273.0_rk)*BB(nloop)*wc(iw)**lp(nloop)
                     enddo
                     wctmp(iwc) = wc(iw)
                  ENDIF

                  IF ((tin(iz) .GE. 210._rk).AND.(tin(iz) .LE. 300._rk)) THEN
                     do nloop = 1,5

                        ytmp(iz,iwc) = ytmp(iz,iwc) &
                             +  AA(nloop)* (wc(iw)**lp(nloop)) &
                             + (tin(iz)-273.0_rk)*BB(nloop)*wc(iw)**lp(nloop) 
                     enddo
                     wctmp(iwc) = wc(iw)
                  ENDIF

                  IF (tin(iz) .GT. 300._rk) THEN
                     do nloop = 1, 5
                        ytmp(iz,iwc) = ytmp(iz,iwc) &
                             +  AA(nloop)* (wc(iw)**lp(nloop)) &
                             + (300.0_rk-273.0_rk)*BB(nloop)*wc(iw)**lp(nloop)
                     enddo
                     wctmp(iwc) = wc(iw)
                  ENDIF
               enddo
               iwc = iwc+ 1

            ENDIF

         enddo

         !     ... Combine cross sections
         do iz = 1, nz
            icnt = 1

            !     ... < 210nm
            do i = 1, n
               IF (xin(i) .LT. 210._rk) THEN
                  ycomb(iz,icnt) = yin(i)
                  wcb  (icnt)    = xin(i)
                  icnt = icnt + 1
               ENDIF
            enddo
            !     ... 210-290 nm
            do i = 1, iwc-1
               ycomb(iz,icnt) = 10**(ytmp(iz,i))
               wcb  (icnt)    = wctmp(i)
               icnt = icnt+1
            enddo
            !     ... >290nm
            do i = 1, n
               IF (xin(i) .GT. 290._rk) THEN
                  ycomb(iz,icnt) = yin(i)
                  wcb  (icnt)    = xin(i)
                  icnt = icnt+1
               ENDIF
            enddo
         enddo
         !----------------------------------------------
         !     ... interpolate
         !----------------------------------------------
         do iz = 1, nz
            n1 = icnt-1
            y1 = ycomb(iz,:)
            x1 = wcb

            call add_pnts_inter2(x1,y1,yg1, kdata, n, &
                 nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)

            if (errflg /= 0) then
               return
            end if

            ytd(iz,:) = yg1(:)

         enddo

         !----------------------------------------------
         !     ...quantum yield assumed to be unity
         !----------------------------------------------
         qy = 1._rk

         do iw = 1, nw-1
            do iz = 1, nz
               sq(iz,iw) = qy * ytd(iz,iw)
            enddo
         enddo
      endif

    contains

      subroutine readit
        x1 = xnan
        y1 = xnan
        n = xsqy_tab(j)%filespec%nread(1)
        call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
             errmsg=errmsg, errflg=errflg, &
             skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
             rd_cnt=n, &
             x=xin,y=yin )

      end subroutine readit

    end subroutine XSQY_CH2BR2

    SUBROUTINE XSQY_MMILLS(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!---------------------------------------------------------------------------!
!  PARAMETERS:                                                              !
!  NW     - INTEGER, number of specified intervals + 1 in working        (I)!
!           wavelength grid                                                 !
!  WL     - REAL, vector of lower limits of wavelength intervals in      (I)!
!           working wavelength grid                                         !
!  WC     - REAL, vector of center points of wavelength intervals in     (I)!
!           working wavelength grid                                         !
!  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)!
!  TLEV   - REAL, temperature (K) at each specified altitude level       (I)!
!  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)!
!  J      - INTEGER, counter for number of weighting functions defined  (IO)!
!  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)!
!           photolysis reaction defined, at each defined wavelength and     !
!           at each defined altitude level                                  !
!  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)!
!           defined                                                         !
!---------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      INTEGER, intent(in) :: nw
      INTEGER, intent(in) :: nz
      INTEGER, intent(inout) :: j
      REAL(rk), intent(in)    :: wl(:), wc(:)
      REAL(rk), intent(in)    :: tlev(:)
      REAL(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

!---------------------------------------------------------------------------!
!     ... local                                                             !
!---------------------------------------------------------------------------!
      INTEGER, parameter :: kdata = 300
      integer, parameter :: jdata = 200
      INTEGER :: n, iw
      REAL(rk) :: x_min(kdata), x_max(kdata), x(kdata), y(kdata)
      REAL(rk) :: yg(kw)
      REAL(rk) :: qy
      real(rk), save :: xsq(kw,jdata)

      if (present(sq)) then
         sq = xnan
      end if
      if( initialize ) then
         xsq(:,j) = xnan
         CALL readit
         xsq(:nw-1,j) = xsqy_tab(j)%qyld*yg(:nw-1)
      else

         sq(1:nw-1,1) = xsq(1:nw-1,j)

      endif

      contains

      subroutine readit
        x_max = xnan
        x_min = xnan
        x = xnan
        y = xnan
        n = xsqy_tab(j)%filespec%nread(1)
        call base_read( filespec=xsqy_tab(j)%filespec%filename(1), &
             errmsg=errmsg, errflg=errflg, &
             skip_cnt=xsqy_tab(j)%filespec%nskip(1), &
             rd_cnt=n, &
             x=x_min, y=x_max, y1=y )
        x(1:n) = 0.5_rk*(x_min(1:n)+x_max(1:n))        
        call add_pnts_inter2(x,y,yg,kdata,n, &
             nw,wl,xsqy_tab(j)%equation,deltax,(/0._rk,0._rk/), errmsg, errflg)
      end subroutine readit

    end subroutine XSQY_MMILLS
    
    subroutine qy_o2(nw,wl,wc,nz,tlev,airden,j, errmsg, errflg, sq )
!---------------------------------------------------------------------------!
!  PARAMETERS:                                                              !
!  NW     - INTEGER, number of specified intervals + 1 in working        (I)!
!           wavelength grid                                                 !
!  WL     - REAL, vector of lower limits of wavelength intervals in      (I)!
!           working wavelength grid                                         !
!  WC     - REAL, vector of center points of wavelength intervals in     (I)!
!           working wavelength grid                                         !
!  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)!
!  TLEV   - REAL, temperature (K) at each specified altitude level       (I)!
!  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)!
!  J      - INTEGER, counter for number of weighting functions defined  (IO)!
!  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)!
!           photolysis reaction defined, at each defined wavelength and     !
!           at each defined altitude level                                  !
!  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)!
!           defined                                                         !
!---------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
!     ... args                                                                !
!-----------------------------------------------------------------------------!
      integer, intent(in) :: nw
      integer, intent(in) :: nz
      integer, intent(inout) :: j
      real(rk), intent(in)    :: wl(:), wc(:)
      real(rk), intent(in)    :: tlev(:)
      real(rk), intent(in)    :: airden(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      real(rk), optional, intent(out) :: sq(:,:)

      integer, save :: srb_ndx
      integer :: iw
      logical, save :: is_initialized = .false.
      
      if (present(sq)) then
         sq = xnan
      end if
!-----------------------------------------------------------------------------
!     ... Shortward of 174.65 the product is O2 + hv => O(3P) + O(1D)
!     ... Longward  of 174.65 the product is O2 + hv => O(3P) + O(3P)
!-----------------------------------------------------------------------------
      if( initialize ) then
         if( .not. is_initialized ) then
            find_srb: do iw = 1, nw-1
               if (wc(iw) > 174.65_rk) then
                  srb_ndx = iw
                  exit find_srb
               end if
            end do find_srb
         end if
      else
!-----------------------------------------------------------------------------
!     ... O2 + hv -> O(3P) + O(1D) at lyman alpha has a qy = 0.53
!         Lacoursiere et al., J. Chem. Phys. 110., 1949-1958, 1999.
!-----------------------------------------------------------------------------
         sq(:nw-1,1) = 0._rk
         
         if (xsqy_tab(j)%channel==1) then
            sq(:srb_ndx-1,1) = 1._rk
            sq(la_ndx,1) = 0.53_rk
         else
            sq(srb_ndx:nw-1,1) = 1._rk
            sq(la_ndx,1) = 0.47_rk
         end if
      end if
      
    end subroutine qy_o2

      SUBROUTINE add_pnts_inter2(xin,yin,yout,kdata,n,nw,wl,jlabel,deltax,yends, errmsg, errflg)

      integer, intent(in) :: kdata
      integer, intent(in) :: n
      integer, intent(in) :: nw
      real(rk), intent(in)    :: deltax
      real(rk), intent(in)    :: wl(nw)
      real(rk), intent(in)    :: xin(kdata)
      real(rk), intent(in)    :: yin(kdata)
      real(rk), intent(in)    :: yends(2)
      real(rk), intent(inout) :: yout(nw-1)
      character(len=*), intent(in) :: jlabel
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: m
      real(rk)    :: xwrk(kdata), ywrk(kdata)

      errmsg = ' '
      errflg = 0

      m = n 
      xwrk(1:n) = xin(1:n)
      ywrk(1:n) = yin(1:n)
      CALL addpnt(xwrk,ywrk,kdata,m,xin(1)*(1._rk-deltax),yends(1),errmsg, errflg)
      CALL addpnt(xwrk,ywrk,kdata,m,              0._rk,yends(1),errmsg, errflg)
      CALL addpnt(xwrk,ywrk,kdata,m,xin(n)*(1._rk+deltax),yends(2),errmsg, errflg)
      CALL addpnt(xwrk,ywrk,kdata,m,          1.e+38_rk,yends(2),errmsg, errflg)

      CALL inter2(nw,wl,yout,m,xwrk,ywrk,errmsg, errflg)

      END SUBROUTINE add_pnts_inter2

      SUBROUTINE base_read( filespec, errmsg, errflg, skip_cnt, rd_cnt,x, y, y1, y2, y3, y4, y5 )

      integer, optional, intent(in) :: skip_cnt
      integer, intent(inout)        :: rd_cnt
      real(rk), intent(inout)           :: x(:), y(:)
      real(rk), optional, intent(inout) :: y1(:), y2(:), y3(:)
      real(rk), optional, intent(inout) :: y4(:), y5(:)
      character(len=*), intent(in)  :: filespec
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, idum
      integer :: y_to_rd
      integer :: ios

      errmsg = ' '
      errflg = 0

      y_to_rd = 1
      if( present(y5) ) y_to_rd = y_to_rd + 1
      if( present(y4) ) y_to_rd = y_to_rd + 1
      if( present(y3) ) y_to_rd = y_to_rd + 1
      if( present(y2) ) y_to_rd = y_to_rd + 1
      if( present(y1) ) y_to_rd = y_to_rd + 1

      OPEN(UNIT=kin,FILE=trim(filespec),STATUS='old',IOSTAT=ios)
      IF( ios /= 0 ) then
         write(errmsg,'(''base_read: failed to open '',a)') trim(filespec)
         errflg = ios
         return
      ENDIF

      if( present(skip_cnt) ) then
        DO i = 1, skip_cnt
          READ(kin,*,IOSTAT=ios)
          IF( ios /= 0 ) exit
        END DO
      else
        READ(kin,*,IOSTAT=ios) idum,rd_cnt
        IF( ios == 0 ) then
          DO i = 1, idum-2
            READ(kin,*,IOSTAT=ios)
            IF( ios /= 0 ) exit
          ENDDO
        ENDIF
      endif

      IF( ios /= 0 ) then
         write(errmsg,'(''base_read: failed to read '',a)') trim(filespec)
         errflg = ios
         return
      ENDIF

      select case( y_to_rd )
        case( 1 )
          DO i = 1, rd_cnt
            READ(kin,*,IOSTAT=ios) x(i), y(i)
            IF( ios /= 0 ) exit
          END DO
        case( 2 )
          DO i = 1, rd_cnt
            READ(kin,*,IOSTAT=ios) x(i), y(i), y1(i)
            IF( ios /= 0 ) exit
          END DO
        case( 3 )
          DO i = 1, rd_cnt
            READ(kin,*,IOSTAT=ios) x(i), y(i), y1(i), y2(i)
            IF( ios /= 0 ) exit
          END DO
        case( 4 )
          DO i = 1, rd_cnt
            READ(kin,*,IOSTAT=ios) x(i), y(i), y1(i), y2(i), y3(i)
            IF( ios /= 0 ) exit
          END DO
        case( 5 )
          DO i = 1, rd_cnt
            READ(kin,*,IOSTAT=ios) x(i), y(i), y1(i), y2(i), y3(i),y4(i)
            IF( ios /= 0 ) exit
          END DO
        case( 6 )
          DO i = 1, rd_cnt
            READ(kin,*,IOSTAT=ios) x(i), y(i), y1(i), y2(i), y3(i),y4(i),y5(i)
            IF( ios /= 0 ) exit
          END DO
      end select

      CLOSE (kin)

      IF( ios /= 0 ) then
         write(errmsg,'(''base_read: failed to read '',a)') trim(filespec)
         errflg = ios
         return
      ENDIF

      END SUBROUTINE base_read

      SUBROUTINE fo3qy2(nz, w, t, qyld)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
! function to calculate the quantum yield O3 + hv -> O(1D) + O2,             =*
! according to:                                                             
! Matsumi, Y., F. J. Comes, G. Hancock, A. Hofzumanhays, A. J. Hynes,
! M. Kawasaki, and A. R. Ravishankara, QUantum yields for production of O(1D)
! in the ultraviolet photolysis of ozone:  Recommendation based on evaluation
! of laboratory data, J. Geophys. Res., 107, 10.1029/2001JD000510, 2002.
!-----------------------------------------------------------------------------*

      INTEGER, intent(in) :: nz
      REAL(rk), intent(in)    :: w
      REAL(rk), intent(in)    :: t(:)
      REAL(rk), intent(inout) :: qyld(:)

      REAL(rk), parameter :: A(3)  = (/ 0.8036_rk, 8.9061_rk, 0.1192_rk/)
      REAL(rk), parameter :: X(3)  = (/ 304.225_rk, 314.957_rk, 310.737_rk/)
      REAL(rk), parameter :: om(3) = (/ 5.576_rk, 6.601_rk, 2.187_rk/)

      REAL(rk), parameter :: q1 = 1._rk

      REAL(rk) :: kt(nz)
      REAL(rk) :: q2(nz), qdiv(nz)

      
      kT(1:nz) = 0.695_rk * t(1:nz)
      q2(1:nz) = exp(-825.518_rk/kT(1:nz))

      kT(1:nz) = t(1:nz)/300._rk
      qdiv(1:nz) = 1/(q1 + q2(1:nz))
      
      IF(w .LE. 305._rk) THEN
        qyld(1:nz) = 0.90_rk
      ELSEIF(w .GT. 305._rk .AND. w .LE. 328._rk) THEN
        qyld(1:nz) = 0.0765_rk + a(1)*q1*qdiv(1:nz)*EXP(-((x(1) - w)/om(1))**4) &
                   + kT(1:nz)*(a(2)*kT(1:nz)*q2*qdiv(1:nz)*EXP(-((x(2) - w)/om(2))**2) &
                               + a(3)*sqrt(kT(1:nz))*EXP(-((x(3) - w)/om(3))**2))
      ELSEIF(w .GT. 328._rk .AND. w .LE. 340._rk) THEN
         qyld(1:nz) = 0.08_rk
      ELSEIF(w .GT. 340._rk) THEN
         qyld(1:nz) = 0._rk
      ENDIF

      END SUBROUTINE fo3qy2

      SUBROUTINE qyacet(nz, w, T, M, fac)
! This file contains subroutines used for calculation of quantum yields for 
! various photoreactions:
!     qyacet - q.y. for acetone, based on Blitz et al. (2004)

! Compute acetone quantum yields according to the parameterization of:
! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield 
!       (2004), Pressure and temperature-dependent quantum yields for the 
!       photodissociation of acetone between 279 and 327.5 nm, Geophys. 
!       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.

      IMPLICIT NONE

! input:
! w = wavelength, nm
! T = temperature, K
! m = air number density, molec. cm-3

      INTEGER, intent(in) :: nz
      REAL(rk), intent(in)    :: w
      REAL(rk), intent(in)    :: T(:), M(:)
      REAL(rk), intent(inout) :: fac(:)

! internal:

      REAL(rk) :: wfac
      REAL(rk) :: a0(nz), a1(nz), a2(nz), a3(nz), a4(nz)
      REAL(rk) :: b0(nz), b1(nz), b2(nz), b3(nz), b4(nz)
      REAL(rk) :: c3(nz)
      REAL(rk) :: cA0(nz), cA1(nz), cA2(nz), cA3(nz), cA4(nz)
      real(rk) :: dumexp(nz)

! fac = quantum yield for product CH3CO (acetyl radical)

      REAL(rk) :: fco(nz)
      REAL(rk) :: tfac(nz)

!** set out-of-range values:
! use low pressure limits for shorter wavelengths
! set to zero beyound 327.5

      IF(w .LT. 279._rk) THEN
        fac(1:nz) = 0.95_rk
      ELSEIF(w .GT. 327._rk) THEN
        fac(1:nz) = 0._rk
      ELSE
        wfac = 1.e7_rk/w
!** CO (carbon monoxide) quantum yields:
        tfac(1:nz) = t(1:nz)/295._rk
        a0(1:nz) = 0.350_rk * tfac(1:nz)**(-1.28_rk)
        b0(1:nz) = 0.068_rk * tfac(1:nz)**(-2.65_rk)
!*SM: prevent exponent overflow in rare cases:

        dumexp(1:nz) = b0(1:nz)*(w - 248._rk)
        where( dumexp(1:nz) > 80._rk )
          cA0(1:nz) = 5.e34_rk
        elsewhere
          cA0(1:nz) = exp(dumexp(1:nz)) * a0(1:nz) / (1._rk - a0(1:nz))
        endwhere

        fco(1:nz) = 1._rk / (1._rk + cA0(1:nz))

!** CH3CO (acetyl radical) quantum yields:

        IF(w >= 279._rk .AND. w < 302._rk) THEN
          a1(1:nz) = 1.600E-19_rk * tfac(1:nz)**(-2.38_rk)
          b1(1:nz) = 0.55E-3_rk   * tfac(1:nz)**(-3.19_rk)
          cA1(1:nz) = a1(1:nz) * EXP(-b1(1:nz)*(wfac - 33113._rk))
          fac(1:nz) = (1._rk - fco(1:nz)) / (1._rk + cA1(1:nz) * M(1:nz))
        ELSEIF(w >= 302._rk .AND. w <= 327._rk) THEN
         a2(1:nz) = 1.62E-17_rk * tfac(1:nz)**(-10.03_rk)
         b2(1:nz) = 1.79E-3_rk  * tfac(1:nz)**(-1.364_rk)
         cA2(1:nz) = a2(1:nz) * EXP(-b2(1:nz)*(wfac - 30488._rk))

         a3(1:nz) = 26.29_rk   * tfac(1:nz)**(-6.59_rk)
         b3(1:nz) = 5.72E-7_rk * tfac(1:nz)**(-2.93_rk)
         c3(1:nz) = 30006._rk  * tfac(1:nz)**(-0.064_rk)
         ca3(1:nz) = a3(1:nz) * EXP(-b3(1:nz)*((1.e7_rk/w) - c3(1:nz))**2)

         a4(1:nz) = 1.67E-15_rk * tfac(1:nz)**(-7.25_rk)
         b4(1:nz) = 2.08E-3_rk  * tfac(1:nz)**(-1.16_rk)
         cA4(1:nz) = a4(1:nz) * EXP(-b4(1:nz)*(wfac - 30488._rk))

         fac(1:nz) = (1._rk - fco(1:nz)) * (1._rk + cA3(1:nz) + cA4(1:nz) * M(1:nz)) &
                     / ((1._rk + cA3(1:nz) + cA2(1:nz) * M(1:nz)) * (1._rk + cA4(1:nz) * M(1:nz)))
        ENDIF
      ENDIF

      END SUBROUTINE qyacet

      SUBROUTINE diagnostics

      integer :: m, n, n1

      open( unit=44,file='TUV.diags')

      write(44,*) 'Photolysis diags'
      write(44,*) ' '
      write(44,'(i3,'' Total photorates'')') npht_tab
      write(44,*) ' '
      do m = 2,npht_tab
        write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
      enddo
      write(44,*) ' '
      write(44,'(''Wrf labels'')')
      write(44,*) ' '
      do m = 2,npht_tab
        write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%rxn_name)
      enddo

      write(44,*) ' '
      write(44,'(i3,'' Photorate(s) with no p,temp dependence'')') &
              count(xsqy_tab(2:npht_tab)%tpflag == 0)
      write(44,*) ' '
      do m = 2,npht_tab
        if( xsqy_tab(m)%tpflag == 0 ) then
          write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
        endif
      enddo

      write(44,*) ' '
      write(44,'(i3,'' Photorate(s) with temp dependence'')') &
              count(xsqy_tab(2:npht_tab)%tpflag == 1)
      write(44,*) ' '
      do m = 2,npht_tab
        if( xsqy_tab(m)%tpflag == 1 ) then
          write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
        endif
      enddo

      write(44,*) ' '
      write(44,'(i3,'' Photorate(s) with press dependence'')') &
              count(xsqy_tab(2:npht_tab)%tpflag == 2)
      write(44,*) ' '
      do m = 2,npht_tab
        if( xsqy_tab(m)%tpflag == 2 ) then
          write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
        endif
      enddo

      write(44,*) ' '
      write(44,'(i3,'' Photorate(s) with temp,press dependence'')') &
              count(xsqy_tab(2:npht_tab)%tpflag == 3)
      write(44,*) ' '
      do m = 2,npht_tab
        if( xsqy_tab(m)%tpflag == 3 ) then
          write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
        endif
      enddo

      write(44,*) ' '
      write(44,'(i3,'' Photorate(s) with second channel'')') &
              count(xsqy_tab(2:npht_tab)%channel == 2)
      write(44,*) ' '
      do m = 2,npht_tab
        if( xsqy_tab(m)%channel == 2 ) then
          write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
        endif
      enddo

      write(44,*) ' '
      write(44,'(i3,'' Photorate(s) with third channel'')') &
              count(xsqy_tab(2:npht_tab)%channel == 3)
      write(44,*) ' '
      do m = 2,npht_tab
        if( xsqy_tab(m)%channel == 3 ) then
          write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
        endif
      enddo

      write(44,*) ' '
      write(44,'(i3,'' Photorate(s) with multiple input files'')') &
              count(xsqy_tab(2:npht_tab)%filespec%nfiles > 1)
      write(44,*) ' '
      do m = 2,npht_tab
        if( xsqy_tab(m)%filespec%nfiles > 1 ) then
          write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
        endif
      enddo

      write(44,*) ' '
      write(44,'('' Photorate(s) with skip == -1'')')
      write(44,*) ' '
      do m = 2,npht_tab
        n = xsqy_tab(m)%filespec%nfiles
        do n1 = 1,n
          if( xsqy_tab(m)%filespec%nskip(n1)  == -1 ) then
            write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
          endif
        enddo
      enddo

      write(44,*) ' '
      write(44,'('' Photorate(s) with skip >= 0'')')
      write(44,*) ' '
      do m = 2,npht_tab
        n = xsqy_tab(m)%filespec%nfiles
        do n1 = 1,n
          if( xsqy_tab(m)%filespec%nskip(n1) >= 0 .and. &
              xsqy_tab(m)%filespec%filename(n1) /= ' ' ) then
            write(44,'(i3,2x,a)') m,trim(xsqy_tab(m)%equation)
          endif
        enddo
      enddo

      write(44,*) ' '
      write(44,'('' Photorate(s) with xfac /= 1.e-20'')')
      write(44,*) ' '
      do m = 2,npht_tab
        n = xsqy_tab(m)%filespec%nfiles
        do n1 = 1,n
          if( xsqy_tab(m)%filespec%xfac(n1) /= 1.e-20_rk ) then
            write(44,'(i3,2x,a,1pg15.7)') &
              m,trim(xsqy_tab(m)%equation),xsqy_tab(m)%filespec%xfac(n1)
          endif
        enddo
      enddo

      write(44,*) ' '
      write(44,'('' Filenames'')')
      write(44,*) ' '
      do m = 2,npht_tab
        n = xsqy_tab(m)%filespec%nfiles
        do n1 = 1,n
          if( xsqy_tab(m)%filespec%filename(n1) /= ' ' ) then
            write(44,'(i3,2x,a,3x,i4,3x,i4)') &
               m,trim(xsqy_tab(m)%filespec%filename(n1)), &
               xsqy_tab(m)%filespec%nskip(n1), &
               xsqy_tab(m)%filespec%nread(n1)
          endif
        enddo
      enddo

      close( 44 )

      END SUBROUTINE diagnostics

      INTEGER FUNCTION get_xsqy_tab_ndx( jlabel,rxn_name )

      character(len=*), optional, intent(in) :: jlabel
      character(len=*), optional, intent(in) :: rxn_name

      integer :: m

      get_xsqy_tab_ndx = -1

      if( present(jlabel) ) then
        do m = 2,npht_tab
          if( trim(jlabel) == trim(xsqy_tab(m)%equation) ) then
            get_xsqy_tab_ndx = m
            exit
          endif
        enddo
      elseif( present(rxn_name) ) then
        do m = 2,npht_tab
          if( trim(rxn_name) == trim(xsqy_tab(m)%rxn_name) ) then
            get_xsqy_tab_ndx = m
            exit
          endif
        enddo
      endif


      END FUNCTION get_xsqy_tab_ndx

      end module module_rxn
