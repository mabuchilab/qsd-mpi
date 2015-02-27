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
! MODULE MESH
! ===========
!
! This module should defines the mesh, that is, the simulation arena in a 
! real-space representation based scheme. It should contain:
! (i) the necessary parameters that define the grid (i.e. number of points, 
!     grid spacing...),
! (ii) functions to relate the "indexes" of the arrays that define each
!     function to real coordinates,
! (iii) the definition of the Hilbert space, which amounts to defining a
!     dot product between functions defined on the grid, and
! (iv) the needed differential operators, which in this case is only the
!     laplacian.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mesh
  implicit none



!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First of all, we need to define the dimensions of the simulation cell, and
! the grid spacing. The next variables are suitable for the examples that will
! be done later.
!
! We will use a regular rectangular mesh, of equal sides (that is, a square).
! This can be easily extended to more complicated geometries, but let us
! keep it simple for the moment.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8), parameter :: mesh_length = 50.0_8 ! the length L of the cube
  integer, parameter :: n = 81       ! total number of points in each direction
  real(8), parameter :: delta = mesh_length/(N-1)  ! spacing between the points




contains




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next functions assign the real "x" and "y" coordinate to a set of integer
! indexes. The following definitions maps the indexes onto the square
! [-L/2, L/2]. It is better if we assume that N is an odd number, and in this
! way (0, 0) belongs to the mesh, and we have an equal number points in each
! direction.
!
! These definitions are once again just simple-minded ones; one could define
! more elabore, non-regular meshes. Note also that we make "x" depend on
! both ix and iy, which is not necessary for this simple example.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8) function x(ix, iy)
   integer, intent(in) :: ix, iy
   x = (ix-1)*delta - (N/2)*delta
end function

real(8) function y(ix, iy)
   integer, intent(in) :: ix, iy
   y = (iy-1)*delta - (N/2)*delta
end function




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To define the Hilbert space necessary for any Quantum-Mechanical code, we
! need to define a dot product. This defines the norm and the distance.
!
! Note that we have two dot products, one for real functions and one for
! complex functions. If you know some Fortran 90, you can bundle the two
! definitions together by means of an module interface; in this way you do not
! have to worry about how the functions are in the rest of the code.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
real(8) function dotproduct(a, b) result(r)
  real(8), intent(in) :: a(n, n), b(n, n)
  r = sum(a(:,:)*b(:, :))*delta**2
end function dotproduct

complex(8) function zdotproduct(a, b) result(r)
  complex(8), intent(in) :: a(n, n), b(n, n)
  r = sum(conjg(a(:,:))*b(:, :))*delta**2
end function zdotproduct




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE LAPLACIAN.
! =====================
!
! INPUT:  
!   f [real(8), dimension(n, n)] : the function whose Laplacian is calculated
! ---------
! OUTPUT: 
!   lapl [real(8), dimension(n, n)] : the Laplacian of f.
!
! The kinetic operator is a Laplacian in real space, so we need a procedure
! that calculates the Laplacian of a function.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine laplacian(f, lapl)

  implicit none

  real(8), intent(in)  :: f(N, N)    ! function whose laplacian is to be taken
  real(8), intent(out) :: lapl(N, N) ! the laplacian


!!!!!! MISSING CODE 1
  integer :: ix, iy, k

  integer, parameter :: order = 4
  real(8), allocatable :: c(:)

  allocate(c(-order:order))
  c(-order:order) = (/ &  ! The coefficients of the laplacian...
      -1.785714d-3, 2.539683d-2, -0.2d0, 1.6d0,      &
      -2.847222d0,                                   &
      1.6d0, -0.2d0, 2.539683d-2, -1.785714d-3 /)

  do ix = 1, N
    do iy = 1, N
      lapl(ix, iy) = 0.0_8
      do k = -order, order
        if(iy+k>=1 .and. iy+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix, iy+k)
        if(ix+k>=1 .and. ix+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix+k, iy)
      end do
    end do
  end do

  lapl = lapl/delta**2
!!!!!! END OF MISSING CODE

end subroutine laplacian




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE ZLAPLACIAN.
! =====================
!
! INPUT:  
!   f [complex(8), dimension(n, n)] : the function whose Laplacian is calculated
! ---------
! OUTPUT: 
!   lapl [complex(8), dimension(n, n)] : the Laplacian of f.
!
! This is exactly the same that laplacian, but for complex functions. The
! missing code is identical.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zlaplacian(f, lapl)

  implicit none

  complex(8), intent(in)  :: f(N, N)    ! function whose laplacian is to be taken
  complex(8), intent(out) :: lapl(N, N) ! the laplacian

!!!!!! MISSING CODE 1
  integer :: ix, iy, k

  integer, parameter :: order = 4
  real(8), allocatable :: c(:)

  allocate(c(-order:order))
  c(-order:order) = (/ &  ! The coefficients of the laplacian...
      -1.785714d-3, 2.539683d-2, -0.2d0, 1.6d0,      &
      -2.847222d0,                                   &
      1.6d0, -0.2d0, 2.539683d-2, -1.785714d-3 /)

  do ix = 1, N
    do iy = 1, N
      lapl(ix, iy) = 0.0_8
      do k = -order, order
        if(iy+k>=1 .and. iy+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix, iy+k)
        if(ix+k>=1 .and. ix+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix+k, iy)
      end do
    end do
  end do

  lapl = lapl/delta**2
!!!!!! END OF MISSING CODE

end subroutine zlaplacian




end module mesh





