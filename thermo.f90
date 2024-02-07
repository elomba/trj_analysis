module mod_thermo
   use mod_precision
   use mod_common
   use mod_input
   use cudafor
   implicit none

contains

   subroutine thermo_calc(iconf)
      integer, intent(in) :: iconf
      integer :: i

      densty = natms/volumen

      If (keytrj > 0) Then
         ekin = masa(natms)*Dot_product(vel(1:ndim, natms&
              &), vel(1:ndim, natms))
      End If
      Do i = 1, natms - 1
         If (keytrj > 0) Then
            ekin = ekin + masa(i)*Dot_product(vel(1:ndim, i)&
                 &, vel(1:ndim, i))
         End if
      End Do
      ekin = 0.5d0*ekin
      ecaver = ecaver + ekin
      If (keytrj > 0) Then
         Tfact = (natms - nmzero - 1)*ndim
         Tfact = (natms - nmzero - 1)*ndim
         taver = 2*ecaver*(aunit/tunit)**2/(Tfact*Rgas*Iconf)
         temperature = 2*ekin*(aunit/tunit)**2/(Tfact*Rgas)
      End if

   end subroutine thermo_calc

end module mod_thermo
