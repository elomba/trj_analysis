!===============================================================================
! Module: mod_precision
!===============================================================================
! Purpose:
!   Defines precision parameters for floating-point arithmetic throughout
!   the trajectory analysis program. Provides flexible precision control
!   by mapping standard Fortran real kinds to custom type names.
!
! Parameters:
!   myprec  - Main precision used for most calculations (single precision)
!   float   - Single precision floating point (real32)
!   double  - Double precision floating point (real64)
!
! Notes:
!   - Default precision is single (sp/real32) for GPU efficiency
!   - Change myprec to dp for double precision calculations if needed
!   - All modules should use these parameters for consistency
!===============================================================================
module mod_precision
    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
    implicit none
    integer, parameter :: myprec = sp
    integer, parameter :: float = sp
    integer, parameter :: double = dp
contains
end module mod_precision
