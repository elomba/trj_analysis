Module fftw3
   Use, Intrinsic :: iso_c_binding
   Include "fftw3.f03"
   Use mod_precision
   Real(float), Parameter :: pi=3.141592653589793
   Real(float) :: dt,dq

Contains

subroutine fftw1d(fin,fout,w,nin,dt, tmax)
   !
   Implicit None
   integer, intent(IN) :: nin
 !  Integer, parameter :: N=512
   real(float), intent(IN) :: fin(nin), dt, tmax
   real(float), intent(OUT) :: fout(nin), w(nin)
   Complex(float), Dimension(:) :: out, in
   Integer :: i, j, nu
   Integer(kind=8) :: plan
   nu = int(log(real(nin))/log(2.00)+0.540)
   n = 2**(nu+1)
   dq = 2*pi/(n*dt)
   allocate(in(n),out(n))
   in(1:nin) =  cmplx(fin(1:nin),0.0)
   in(nin+1:n) = (0.0,0.0)
   
   Call fftw2(in,out,n,.True.)
   Do i=1,nin
      w(i) = (i-1)*dq
      fout(i) = real(out(i))
    !  Write(388,'(10f15.7)')q, out(i), (out(i)+out(n-i+2))*dr, sqrt(pi/sig12)*exp(-q**2/4/sig12)
   End Do
 end subroutine fftw1d

  Subroutine fftw1(in,out,n,forward)
   Implicit None
   Integer, Intent(IN) :: n
   Complex(float), Dimension(n), Intent(INOUT) :: in, out
   Logical, Intent(in) :: forward
   Real (float) :: dfact
   Integer(float) :: plan
   Integer :: direction, i, j
   If (forward) Then
      dfact = dt
      direction = -1
   Else
      dfact = dq/(2*pi)
      direction = +1
   End If
   !
   ! Correct for trapezoidal rule
   !
   in(1) = 0.5*in(1)
   Call fftw_plan_dft_1d(plan,N,in,out,direction,FFTW_ESTIMATE)
   Call fftw_execute_dft(plan,in,out)
   Call fftw_destroy_plan(plan)
   out(1) = 2*out(1)*dfact
   Forall (i=2:n/2)
      ! FT(k_x,k_y) = 2*(Real(F(k_x,k_y))+Real(F(K_x,-k_y)) and we make use of the symmetry 
      out(i) = (out(i)+out(n-i+2))*dfact
   End Forall
 End Subroutine fftw1

 function window(t,x,n,alpha) return (result)
   implicit none
   integer, intent(in) :: n
   real(float), intent(in) :: t(n), x
   real(flota), intent(out) :: result(n)
   result(1:n) = 0.5*(1-tanh(alpha*(t(1:n)-x)))
 end function window

end module fftw3
