!===============================================================================
! Module: fftw3
!===============================================================================
! Purpose:
!   Provides wrapper routines for FFTW3 library to perform Fast Fourier
!   Transforms. Used for computing dynamic structure factors S(q,ω) by
!   transforming time correlation functions from time to frequency domain.
!
! Key Functionality:
!   - 1D Fourier transforms (time to frequency)
!   - Window function application for finite time series
!   - Automatic zero-padding for FFT efficiency
!   - Trapezoidal rule correction
!
! Main Variables:
!   dt - Time step spacing
!   dq - Frequency spacing (2π/(N*dt))
!   pi - Pi constant
!
! Subroutines:
!   fftw1d()   - High-level 1D FFT with windowing and zero-padding
!   fftw1()    - Low-level FFTW3 wrapper for complex FFT
!   window()   - Tanh window function for smooth truncation
!
! Parameters:
!   fin(:)   - Input time series
!   fout(:)  - Output frequency spectrum
!   w(:)     - Frequency array
!   nin      - Number of input points
!   dtin     - Time step
!   tmax     - Maximum time for windowing
!
! Notes:
!   - Requires FFTW3 library installation
!   - Uses fftw3.f03 include file for FFTW interface
!   - Automatically pads to next power of 2 for efficiency
!   - Window function smoothly truncates data at tmax
!   - Applies trapezoidal rule correction at boundaries
!===============================================================================
Module fftw3
   Use, Intrinsic :: iso_c_binding
   Use mod_precision
   Include "fftw3.f03"
   Real(double), Parameter :: pi=3.14159265358979323846_double
   Real(double) :: dt,dq

Contains


subroutine fftw1d(fin,fout,w,nin,nout,dtin, tmax)
   !
   !
   Implicit None
   integer, intent(IN) :: nin, nout
   real(double), intent(IN) :: fin(nin), dtin, tmax
   real(double), intent(OUT) :: fout(nout), w(nout)
   real(double), allocatable, dimension(:) :: tx
   Complex(double), allocatable, Dimension(:) :: out, in
   Real(double) :: alpha=1.0_double
   Integer :: i, j, nu, n
   Integer(kind=8) :: plan
   ! Determine next power of 2 such that n/2 >= nin
   nu = ceiling(log(real(nin))/log(2.00))
   n = 2**(nu+1)
   ! Set up time and frequency grids
   dt = dtin
   dq = 2*pi/(n*dt)
   allocate(in(n),out(n),tx(n))
   ! Copy input data and zero-pad
   in(1:nin) =  cmplx(fin(1:nin),0.0_double,kind=double)
   in(nin+1:n) = (0.0_double,0.0_double)
   forall (i=1:n) tx(i)=(i-1)*dt
   ! if tmax not defined use maximum time
   if (abs(tmax)<1.0e-6_double) then
      in(1:n) = in(1:n)*window(tx,n*dt,n,alpha)
   else
      in(1:n) = in(1:n)*window(tx,tmax,n,alpha)
   endif

   ! Perform FFT
   Call fftw1(in,out,n,.True.)
   ! Extract positive frequencies and scale
   Do i=1, min(nout, n/2+1)
      w(i) = (i-1)*dq
      fout(i) = real(out(i), kind=double)
   End Do
   if (nout > n/2+1) then
      Do i = n/2+2, nout
         w(i) = (i-1)*dq
         fout(i) = 0.0_double
      End Do
   endif
   deallocate(in,out,tx)
end subroutine fftw1d

Subroutine fftw1(in,out,n,forward)
   !
   ! Use FFTW3 to compute the 1D Fourier cosine transform from t->w
   !
   Implicit None
   Integer, Intent(IN) :: n
   Complex(double), Dimension(n), Intent(INOUT) :: in, out
   Logical, Intent(in) :: forward
   Real (double) :: dfact
   Integer(kind=8) :: plan
   Integer :: direction, i, j
   ! Set transform direction and scaling factor
   If (forward) Then
      dfact = dt
      direction = -1
   Else
      dfact = dq/(2*pi)
      direction = +1
   End If
   !
   ! Correct for trapezoidal rule (first point weighted by 0.5)
   !
   in(1) = 0.5*in(1)
   ! Create FFTW plan and execute
   Call dfftw_plan_dft_1d(plan,N,in,out,direction,FFTW_ESTIMATE)
   Call dfftw_execute_dft(plan,in,out)
   Call dfftw_destroy_plan(plan)
   ! Scale output and use symmetry to extract real part
   out(1) = cmplx(real(out(1))*dfact, 0.0_double, kind=double)
   ! For k=2 to N/2, use symmetry: F(k) = F*(-k) for real input
   Forall (i=2:n/2)
      out(i) = cmplx(real(out(i))*dfact, 0.0_double, kind=double)
   End Forall
   out(n/2+1) = cmplx(real(out(n/2+1))*dfact, 0.0_double, kind=double)
End Subroutine fftw1

function window(t,x,n,alpha) result(w)
   implicit none
   integer, intent(in) :: n
   real(double), intent(in) :: t(n), x, alpha
   real(double) :: w(n)
   w(1:n) = 0.5_double*(1.0_double - tanh(alpha*(t(1:n) - x)))
end function window

end module fftw3
