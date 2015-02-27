!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE POISSON
! ==============
!
! 
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module poisson
use mesh
use cube_function
implicit none




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! There are four public procedures in this module; all the rest is private.
  private
  public :: poisson_init, &
            poisson_solve, &
            poisson_sum,   &
            poisson_fft




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The next set of parameters fixes the way in which the Hartree term is
  ! calculated: performing directly the sum (mode = HARTREE_SUM), or via
  ! the FFT technique (mode = HARTREE_FFT).
  integer, parameter :: HARTREE_SUM = 1, &
                        HARTREE_FFT = 2
  integer, parameter :: mode = HARTREE_FFT




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The next parameters set the type of interaction that defines the system:
  ! Coulombic (interaction = COULOMB), or of Yukawa form (interaction = YUKAWA).
  ! In this latter case, gamma is the parameter that defines the
  ! defines the Yukawa potential.
  integer, parameter :: COULOMB = 1, &
                        YUKAWA  = 2
  integer, parameter :: interaction = COULOMB
  real(8), parameter :: gamma = 2._8




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These variables are only meant for internal use within the module.
  type(dcf) :: fft_cf
  real(8), pointer :: fft_coulb_FS(:,:)
  real(8), parameter :: pi = 3.141592653589793_8

contains




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE POISSON_SOLVE
! ========================
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poisson_solve(rho, v)
  use mesh

  implicit none

  real(8), intent(in)  :: rho(N, N)  ! the density
  real(8), intent(out) :: v(N, N)    ! the solution of the poisson equation

    call poisson_sum(rho, v, N)

end subroutine poisson_solve

subroutine poisson_sum(rho, v, num)
!  use mesh

  implicit none

  integer, intent(in)  :: num        ! grid size
  real(8), intent(in)  :: rho(N, N)  ! the density
  real(8), intent(out) :: v(N, N)    ! the solution of the poisson equation

  integer  :: ix, iy, jx, jy
  real(8)  :: r1(2), r2(2)
  real(8), parameter :: pi = 3.141592653589793d0

!  v = 0.0d0
  write(*, *) 'starting poisson sum', num
  do ix = 1, num
    do iy = 1, num
      r1(1) = x(ix, iy); r1(2) = y(ix, iy)
      do jx = 1, N
        do jy = 1, N
           r2(1) = x(jx, jy); r2(2) = y(jx, jy)

          if(ix == jx .and. iy == jy) then
            v(ix, iy) = v(ix, iy) + 2.0d0*sqrt(pi)*rho(ix, iy)/delta
          else
            v(ix, iy) = v(ix, iy) + rho(jx, jy)/sqrt(sum((r1-r2)**2))
          end if
        end do
      end do
      v(ix, iy) = v(ix, iy)*delta**2
    end do
  end do
  write(*, *) 'Done with Poisson sum'
end subroutine poisson_sum

subroutine poisson_fft(rho, pot, n)
  implicit none
  integer, intent(in)  ::  n
  real(8), intent(in)  ::  rho(n, n)
  real(8), intent(out) ::  pot(n, n)

  integer :: k, ix, iy

  call dcf_alloc_RS(fft_cf)          ! allocate the cube in real space
  call dcf_alloc_FS(fft_cf)          ! allocate the cube in Fourier space
  write(*, *) 'Allocated RS and FS', n
  fft_cf%rs = 0.0_8
  do ix = 1, n
     do iy = 1, n
        fft_cf%rs(fft_cf%n(1)/2+1-n/2+ix,fft_cf%n(2)/2+1-n/2+iy) = rho(ix, iy)
     enddo
     write(*, *) 'rhos assigned ', ix,  fft_cf%rs(fft_cf%n(1)/2+1-n/2+ix,fft_cf%n(2)/2+1-n/2+iy)
  enddo
  
   write(*, *) 'Calling FT', n

  call dcf_RS2FS(fft_cf)             ! Fourier transform

  ! multiply by the FS of the Coulomb interaction
  fft_cf%FS(1:fft_cf%nx,1:fft_cf%n(2)) = fft_cf%FS(1:fft_cf%nx,1:fft_cf%n(2))*&
                                         fft_Coulb_FS(1:fft_cf%nx,1:fft_cf%n(2))


  call dcf_FS2RS(fft_cf)             ! Fourier transform back

  do ix = 1, n
     do iy = 1, n
        pot(ix, iy) = fft_cf%rs(fft_cf%n(1)/2+1-n/2+ix,fft_cf%n(2)/2+1-n/2+iy)
     enddo
  enddo

  call dcf_free_RS(fft_cf)           ! memory is no longer needed
  call dcf_free_FS(fft_cf)

end subroutine poisson_fft

subroutine poisson_init()
  implicit none
    integer :: ix, iy, iz, ixx(2), db(2)
    real(8) :: temp(2), vec
    real(8) :: gpar,gperp,gx,gz,r_c
    real(8) :: DELTA_R = 1.0e-12_8

    ! double the box to perform the fourier transforms
    db(1) = nint((1.0_8+sqrt(2.0_8))*n)
    db(2) = nint((1.0_8+sqrt(2.0_8))*n)
    write(*, *) 'Calling dcf new'
    call dcf_new(db, fft_cf)    ! allocate cube function where we will perform
    write(*, *) 'Calling dcf fft init', delta, n
    call dcf_fft_init(fft_cf)   ! the ffts

    r_c = sqrt(2.0_8)*n*delta

    ! store the fourier transform of the Coulomb interaction
    allocate(fft_Coulb_FS(fft_cf%nx, fft_cf%n(2)))
    fft_Coulb_FS = 0.0_8

    ! Warning: this has to be set....
    temp(1) = 2.0_8*pi/(db(1)*delta)
    temp(2) = 2.0_8*pi/(db(2)*delta)
    write(*, *) 'Calling besselint'
    do ix = 1,1000
          vec = besselint((ix-1)*1.0_8)
          write(73, *) (ix-1)*1.0_8, r_c*vec
    enddo

      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do ix = 1, fft_cf%nx
          ixx(1) = pad_feq(ix, db(1), .true.)
           vec = sqrt(sum((temp(:)*ixx(:))**2))
           select case(interaction)
             case(COULOMB)
!!$               fft_coulb_fs(ix, iy) = r_c*besselint(vec*r_c)
               fft_coulb_fs(ix, iy) = 2.0 * pi * r_c*besselint(vec*r_c)
             case(YUKAWA)
               fft_coulb_fs(ix, iy) = yukawa_fs(vec)
             end select
      end do
      end do
     write(*, *) 'Done poisson initiation'
end subroutine poisson_init


real(8) function yukawa_fs(x) result(y)
  real(8), intent(in) :: x
  y = 2.0_8*pi/(gamma*sqrt(1.0_8+x**2/gamma**2))
end function yukawa_fs

!!$real(8) function besselint(x) result(y)
!!$  real(8), intent(in) :: x
!!$  real(8), external :: bessel
!!$  integer :: k
!!$  real(8) :: z
!!$  if(x < 0.2_8) then
!!$     y = 2*pi - (pi/6.0_8)*x**2
!!$     return
!!$  endif
!!$  y = 0.0_8
!!$  k = 1
!!$  do
!!$    z = bessel(k, x)/x
!!$    y = y + z
!!$    if(abs(z)<1.0e-9_8) exit
!!$    k = k + 2
!!$  enddo
!!$  y = 4*pi*y
!!$end function


  ! ---------------------------------------------------------
  ! F(x) = (1/x) Integrate[ BesselJ[0, r], {0, x, r} ] =
  !      = HypergeometricPFQ[ {1/2}, {1,3/2}, -x*x/r ] =
  !      = (1/x) * 2 * sum_{k=0}^{\infty} BesselJ[k, x]
real(8) function besselint(x) result(y)
  implicit none

  real(8), intent(in) :: x
  integer :: k, nmax
  real(8) :: z, s
  real(8), allocatable :: bess(:)

  real(8), parameter :: large = 1.0e10_8


  if(x < 0.2_8) then
    y = 1.0_8 - (1.0_8/12.0_8)*x**2
    return
  end if

  nmax = 0

  main_loop: do
    nmax = nmax + 100
    if(.not.allocated(bess)) allocate(bess(0:nmax))

    ! We need to do a backwards recursion since otherwise it is unstable.
    bess(0:nmax) = 0.0_8
    bess(nmax) = 0.0_8
    bess(nmax-1) = 1.0_8
    s = bess(nmax)
    do k = nmax - 2, 0, -1
       bess(k) = (2.0_8*(k+1)/x)*bess(k+1) - bess(k+2)
       if(bess(k) > large) then
         bess(k:nmax) = bess(k:nmax) / large
         s = s / large
       end if
       if(mod(k,2).eq.0) s = s + bess(k)
    end do
    s = 2*s - bess(0)
    do k = 0, nmax
      bess(k) = bess(k)/s
    end do

    y = 0.0_8
    k = 2
      ! In the sum, I use the recursion relation to eliminate half of the terms.
    do
      if(k + sqrt(40.0_8*k) > nmax) exit 
         ! Beyond this value, bess(k) could be imprecise.
      if(mod(k-2,4).eq.0) then
        z = 2*k*bess(k)/x**2
        y = y + z
        if(abs(z) < 1.0e-9_8) exit main_loop
      end if
      k = k + 1
    end do

      deallocate(bess)

    end do main_loop

  if(allocated(bess)) deallocate(bess)
  y = 2.0_8*y
end function besselint




end module poisson


