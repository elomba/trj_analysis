!===============================================================================
! Module: g_types
!===============================================================================
! Purpose:
!   Defines precision parameters for the NetCDF reading utilities.
!   Provides consistent floating-point precision across the read_nc
!   and lmptrj2nc tools.
!
! Parameters:
!   double  - Double precision (15 decimal digits)
!   float   - Single precision (6 decimal digits)
!   myprec  - Default precision for calculations (float)
!
! Notes:
!   - Uses selected_real_kind for portable precision specification
!   - Default is single precision for memory efficiency
!===============================================================================
module g_types
  ! Define general types
  integer, parameter :: double = selected_real_kind(15)
  integer, parameter :: float = selected_real_kind(6)
  integer, parameter :: myprec = float
end module g_types
