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

subroutine output(f, filename)
  use mesh
  implicit none

  real(8), intent(in) :: f(N, N)    ! function to output
  character(len=*), intent(in) :: filename ! file to output

  integer, parameter :: iunit = 99
  integer :: ix, iy

  open(iunit, file=trim(filename), status='unknown')
  do ix = 1, N
    do iy = 1, N
      write(iunit, *) x(ix, iy), y(ix, iy), f(ix, iy)
    end do
    write(iunit, '(1x)')
  end do

  close(iunit)
end subroutine output
