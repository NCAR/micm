!> \file solver_var_defs.f90
!!  Contains type definitions for solver variables

module solver_var_defs

 use ODE_solver, only : baseOdeSolver

 implicit none
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! The following definition sets up the variables for use within solver
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 type Solver_type
   class(baseOdeSolver), pointer :: theSolver
 end type Solver_type


contains

end module solver_var_defs
