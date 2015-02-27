subroutine psn(n, rho)
 use mesh
  use fft
  use poisson
  implicit none

  integer, intent(in):: n
  real(8), intent(in) :: rho(n, n)
  real(8), allocatable :: rh(:, :), vhsum(:, :), vhexact(:, :), vhfft(:, :)


  integer   :: ix, iy
  real(8)   :: r2, alpha
  real(8), parameter :: pi = 3.141592653589793_8
  allocate(vhsum(n, n), vhexact(n, n), vhfft(n, n), rh(n, n))

  write(*, *) 'Grid size ', n
   alpha = 3.0
  do ix = 1, n
     do iy = 1, n
        r2 = x(ix, iy)**2 + y(ix, iy)**2
        rh(ix, iy) = exp(-r2/alpha**2)
     enddo
  enddo
  rh(:, :) = rh(:, :)/(alpha**2*pi)

  call fft_all_init()

    write(*, *) 'Calling poisson_init...'
  call poisson_init()
   write(*, *) 'Calling poisson_solve...', rho(2, 2)
  call poisson_sum(rh, vhsum, n)
   write(*, *) 'Done.'
  call output(vhsum, 'vhsum')
   write(*, *) 'Calling poisson_fft...'
  call poisson_fft(rho, vhfft, n)
   write(*, *) 'Done.'
  call output(vhfft, 'vhfft')

   deallocate(vhsum, vhfft, vhexact)

end subroutine psn
