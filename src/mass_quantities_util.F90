!--------------------------------------------------------------------------------
! utility to compute background air mass quantities
!--------------------------------------------------------------------------------
module mass_quantities_util
  use machine, only : rk => kind_phys
  use const_props_mod, only : const_props_type
  implicit none
  
  integer :: o2_ndx, n2_ndx
  real(rk), allocatable :: molar_mass(:)
  
contains

!> \section arg_table_mass_quantities_util_init Argument Table
!! | local_name | standard_name              | long_name                  | units   | rank | type             | kind      | intent | optional |
!! |------------|----------------------------|----------------------------|---------|------|------------------|-----------|--------|----------|
!! | cnst_info  | chemistry_constituent_info | chemistry_constituent_info | DDT     |    1 | const_props_type |           | in     | F        |
!! | errmsg     | ccpp_error_message         | CCPP error message         | none    |    0 | character        | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag            | CCPP error flag            | flag    |    0 | integer          |           | out    | F        |
!!
  subroutine mass_quantities_util_init(cnst_info, errmsg, errflg)
    type(const_props_type), intent(in)  :: cnst_info(:)
    character(len=512),     intent(out) :: errmsg
    integer,                intent(out) :: errflg

    integer :: ncnst, i

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    o2_ndx =-1
    n2_ndx =-1

    ncnst = size(cnst_info)
    allocate(molar_mass(ncnst))

    do i = 1,ncnst
       molar_mass(i) = cnst_info(i)%get_wght()
       if ( cnst_info(i)%get_name() == 'O2' ) o2_ndx = i
       if ( cnst_info(i)%get_name() == 'N2' ) n2_ndx = i
    end do

  end subroutine mass_quantities_util_init

!> \section arg_table_mass_quantities_util_run Argument Table
!! | local_name | standard_name             | long_name                 | units         | rank | type      | kind      | intent | optional |
!! |------------|---------------------------|---------------------------|---------------|------|-----------|-----------|--------|----------|
!! | press      | pressure                  | ambient pressure          | Pa            | 0    | real      | kind_phys | in     | F        |
!! | temp       | temperature               | ambient temperature       | K             | 0    | real      | kind_phys | in     | F        |
!! | vmr        | concentration             | species concentration     | mole/mole     | 1    | real      | kind_phys | in     | F        |
!! | density    | total_number_density      | total number density      | molecules/cm3 | 0    | real      | kind_phys | out    | F        |
!! | mbar       | mean_molec_mass           | mean molecular mass       | g/mole        | 0    | real      | kind_phys | out    | F        |
!! | errmsg     | ccpp_error_message        | CCPP error message        | none          | 0    | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag           | CCPP error flag           | flag          | 0    | integer   |           | out    | F        |
!!
  subroutine mass_quantities_util_run( press, temp, vmr, density, mbar, errmsg, errflg )

    real(rk), intent(in)            :: press
    real(rk), intent(in)            :: temp
    real(rk), intent(in)            :: vmr(:)
    real(rk), intent(out)           :: density
    real(rk), intent(out)           :: mbar
    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errflg

    real(rk) :: n2_vmr
    real(rk), parameter :: molar_mass_n2 = 28.0134_rk ! g/mole
    real(rk), parameter :: kboltz= 1.38064852e-16_rk ! boltzmann constant (erg/K)
   
    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    density = 10._rk*press/(kboltz*temp)

    if (o2_ndx>0) then
       mbar = sum( vmr(:)*molar_mass(:) )
       if (n2_ndx<-1) then
          n2_vmr = 1._rk - sum(vmr(:))
          mbar = mbar + n2_vmr*molar_mass_n2
       endif
    else
       mbar = 28.966_rk ! set to constant if the major species VMRs are unknown
    endif
    
  end subroutine mass_quantities_util_run

!> \section arg_table_mass_quantities_util_finalize Argument Table
!! | local_name | standard_name         | long_name                      | units   | rank | type      | kind      | intent | optional |
!! |------------|-----------------------|--------------------------------|---------|------|-----------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message    | CCPP error message             | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag       | CCPP error flag                | flag    |    0 | integer   |           | out    | F        |
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
