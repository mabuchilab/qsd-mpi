!! Copyright (C) 2002-2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
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
!! 02111-1307, USA.

!!!/*#include "global.h"*/

module cube_function
  use fft

  implicit none

  type dcf
    integer :: n(2)
    real(8), pointer :: RS(:,:)
    complex(8), pointer :: FS(:,:)
    integer :: nx ! = n(1)/2 + 1, first dimension of the FS array
    type(fft_type), pointer :: fft
  end type dcf

contains

subroutine dcf_new(n, cf)
!!$  integer, intent(in) :: n(3)
  integer, intent(in) :: n(2)
  type(dcf), intent(out) :: cf
  nullify(cf%RS)
  nullify(cf%FS)
  cf%n = n
  nullify(cf%fft)
end subroutine dcf_new

subroutine dcf_alloc_RS(cf)
  type(dcf), intent(inout) :: cf
  allocate(cf%RS(cf%n(1), cf%n(2)))
end subroutine dcf_alloc_RS

subroutine dcf_free_RS(cf)
  type(dcf), intent(inout) :: cf
  deallocate(cf%RS)
  nullify(cf%RS)
end subroutine dcf_free_RS

subroutine dcf_alloc_FS(cf)
  type(dcf), intent(inout) :: cf
  allocate(cf%FS(cf%nx, cf%n(2)))
end subroutine dcf_alloc_FS

subroutine dcf_free_FS(cf)
  type(dcf), intent(inout) :: cf
  deallocate(cf%FS)
  nullify(cf%FS)
end subroutine dcf_free_FS

subroutine dcf_free(cf)
  type(dcf), intent(inout) :: cf

  if(associated(cf%RS)) then
    deallocate(cf%RS)
    nullify(cf%RS)
  end if

  if(associated(cf%FS)) then
    deallocate(cf%FS)
    nullify(cf%FS)
  end if

  if(associated(cf%fft)) then
    call fft_end(cf%fft)
    deallocate(cf%fft)
    nullify(cf%fft)
  end if

end subroutine dcf_free

!!$!!! initializes the ffts. As the dimension of the fft may be adjusted, this
!!$!!! routine has to be called before allocating anything
subroutine dcf_fft_init(cf)
  type(dcf), intent(inout)  :: cf
  allocate(cf%fft)
  call fft_init(cf%n, fft_real, cf%fft)
  cf%nx = cf%n(1)/2 + 1
end subroutine dcf_fft_init

!!! The next routines convert the function between real space and fourier space
!!! Note that the dimensions of the function in FS are different whether f 
!!! is real or complex, because the FFT representation is different (FFTW scheme).
subroutine dcf_RS2FS(cf)
  type(dcf), intent(inout)  :: cf

!!$   ASSERT(associated(cf%RS))
!!
  if(.not.associated(cf%FS)) call dcf_alloc_FS(cf)

  call dfft_forward(cf%fft, cf%RS, cf%FS)

end subroutine dcf_RS2FS

subroutine dcf_FS2RS(cf)
  type(dcf), intent(inout)  :: cf

!!$  ASSERT(associated(cf%FS))
  if(.not.associated(cf%RS)) call dcf_alloc_RS(cf)

  call dfft_backward(cf%fft, cf%FS, cf%RS)

end subroutine dcf_FS2RS


!!$#include "undef.F90"
!!$#include "complex.F90"
!!$#include "cf_inc.F90"

end module cube_function
  

