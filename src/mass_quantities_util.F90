!--------------------------------------------------------------------------------
! utility to compute background air mass quantities
!--------------------------------------------------------------------------------
module mass_quantities_util
!  USE ccpp_kinds, ONLY: rk => kind_phys
  USE ccpp_kinds, ONLY: kind_phys
  use const_props_mod, only : const_props_type
  implicit none
  
  integer :: o2_ndx, n2_ndx
  real(kind_phys), allocatable :: molar_mass(:)
  
contains

!> \section arg_table_mass_quantities_util_init Argument Table
!! \htmlinclude mass_quantities_util_init.html
!!
  subroutine mass_quantities_util_init(cnst_info, ncnst, errmsg, errflg)
    type(const_props_type), intent(in)  :: cnst_info(:)
    integer,                intent(in)  :: ncnst
    character(len=512),     intent(out) :: errmsg
    integer,                intent(out) :: errflg

    integer :: i

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    o2_ndx =-1
    n2_ndx =-1

    allocate(molar_mass(ncnst))

    do i = 1,ncnst
       molar_mass(i) = cnst_info(i)%get_wght()
       if ( cnst_info(i)%get_name() == 'O2' ) o2_ndx = i
       if ( cnst_info(i)%get_name() == 'N2' ) n2_ndx = i
    end do

  end subroutine mass_quantities_util_init

!> \section arg_table_mass_quantities_util_run Argument Table
!! \htmlinclude mass_quantities_util_run.html
!!
  subroutine mass_quantities_util_run( press, temperature, vmr, density, mbar, errmsg, errflg )

    real(kind_phys), intent(in)            :: press
    real(kind_phys), intent(in)            :: temperature
    real(kind_phys), intent(inout)         :: vmr(:)
    real(kind_phys), intent(out)           :: density
    real(kind_phys), intent(out)           :: mbar
    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errflg

    real(kind_phys) :: n2_vmr
    real(kind_phys), parameter :: molar_mass_n2 = 28.0134_kind_phys ! g/mole
    real(kind_phys), parameter :: kboltz= 1.38064852e-16_kind_phys ! boltzmann constant (erg/K)
   
    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    density = 10._kind_phys*press/(kboltz*temperature)

    if (o2_ndx>0) then
       mbar = sum( vmr(:)*molar_mass(:) )
       if (n2_ndx<-1) then
          n2_vmr = 1._kind_phys - sum(vmr(:))
          mbar = mbar + n2_vmr*molar_mass_n2
       endif
    else
       mbar = 28.966_kind_phys ! set to constant if the major species VMRs are unknown
    endif
    
  end subroutine mass_quantities_util_run

!> \section arg_table_mass_quantities_util_finalize Argument Table
!! \htmlinclude mass_quantities_util_finalize.html
!!
  subroutine mass_quantities_util_finalize( errmsg, errflg )

    !--- arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    deallocate(molar_mass)
    
  end subroutine mass_quantities_util_finalize
  
end module mass_quantities_util
