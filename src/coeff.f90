!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAM COEFF
! =============
!
! It outputs the coefficients of the nine-point formula for the second
! derivative in a regular one-dimensional grid (a Laplacian in two or three
! dimensions may be easily worked out by summing the terms).
!
! By changing the variables "m" (derivative order) and "n" (order of the
! approximation, e.g 4 for a 2*4+1=9 points formula), you may obtain the
! coefficients for other cases.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fcoeff(n, c)

  integer, intent(in):: n
  real(8), intent(in) :: c(2*n+1)

  integer :: m, res
!  real(8), allocatable :: c(:), x(:)
   real(8), allocatable :: x(:)
  real(8) :: delta

  delta = 1.0

  m = 2
!!  n = 4

!  allocate(c(2*n+1), x(2*n))
   allocate(x(2*n))

  do i = 1, n
    x(i) = i*delta
    x(n+i) = -1*i*delta
  end do
!  x(1) = delta
!  x(2) = 2*delta
!  x(3) = 3*delta
!  x(4) = 4*delta
!  x(5) = -delta
!  x(6) = -2*delta
!  x(7) = -3*delta
!  x(8) = -4*delta

  call coefficients(m, n, x, c)

  ! Note that the indexes are changed... just notation.
  write(*, '(a7,e16.7)')'c(0)     = ', c(2*n+1)
  write(*, *)
  write(*,'(a,8e16.7)') 'c( 1: n) = ', c(2:n+1)
  write(*, *)
  write(*,'(a,8e16.7)') 'c(-1:-n) = ', c(n+2:2*n+1)

!  open (unit=1, file = 'coeff.dat', status = 'new')
!   write(1, '(a7,e16.7)')'c(0)     = ', c(2*n+1)
!   write(1, *)
!   write(1,'(a,8e16.7)') 'c( 1: n) = ', c(2:n+1)
!   write(1, *)
!   write(1,'(a,8e16.7)') 'c(-1:-n) = ', c(n+2:2*n+1)
end subroutine fcoeff




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION COEFF
! ==============
!
! INPUT:
!   m [integer] : the order of the derivatives representation.
!   n [integer] : The number of points given is 2*n
!   x [real(8), dimension(2*n)] : positions of the points. The "problem" point position is not given,
!     and assumed to be zero.
! ---------
! OUTPUT:
!   c [real(8), dimension(2*n+1)] : the coefficients of the points. The first one corresponds to the
!     the coefficient at the problem points (which is always minus the sum of all the others), whereas
!     the rest are ordered in the same manner that were given in array x.
!   coeff [integer] : error code. It is the error code of the LAPACK subroutine dgels
!
! Calculates the coefficients for the representation of the m-th order
! derivative of a function at a given point, given that we will have access
! to the values of this function at 2*n points around it (besides the value
! of this function at the problem point). Typically this means n points to the
! left, and n points to the right, but that is not mandatory.
!
! NOTES:
! ------
! It requires BLAS dgesv subroutine.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!integer function coefficients(m, n, x, c) result(ierr)
 subroutine coefficients(m, n, x, c)
  integer, intent(in) :: m, n
  real(8), intent(in) :: x(2*n)
  real(8), intent(out) :: c(2*n+1)

  integer :: i, j, k, lwork, info, a1, a2, a3, a4
  real(8), allocatable :: a(:, :), e(:), work(:)

  allocate(a(2*n, 2*n), e(2*n))

  do i = 1, 2*n
     do j = 1, 2*n
        a(i, j) = x(j)**i
     enddo
  enddo

  k = 1
  e = 0.0
  do i = 1, 2*n
     k = k*i
     if(m==i) then
       e(i) = k
       exit
     endif
  enddo

  a1 = 2*n
  a2 = 2*n
  a3 = 2*n
  a4 = 2*n
  lwork = -1
  allocate(work(1))
   print *, ' a4 ', a4, a(5,5)
  call dgels('n', a1, a2, 1, a, a3, e, a4, work, lwork, info)
  lwork = work(1)
  deallocate(work); allocate(work(lwork))
  print *, 'first dgels call over'
  call dgels('n', a1, a2, 1, a, a3, e, a4, work, lwork, info)

  c(1) = - sum(e(1:2*n))
 do j = 1, 2*n
     c(j+1) = e(j)
  enddo

  deallocate(work, a, e)
  ierr = info
!!end function coefficients
end subroutine coefficients
