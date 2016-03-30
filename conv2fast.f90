! fast convolution routine using fortran 90

module conv2fast
  
  implicit none
  include "omp_lib.h"
  
  real, dimension(4100, 4100) :: iout

contains

  ! set number of threads
  subroutine set_num_threads(n)

    implicit none

    !f2py threadsafe
    !f2py intent(in) n

    INTEGER :: n
    CALL OMP_SET_NUM_THREADS(n)

  end subroutine set_num_threads

  ! convolution
  subroutine conv (ninx, niny, Nifx, Nify, nf, kernel, iin)

    implicit none

    !f2py threadsafe
    !f2py intent(in) ninx, niny
    !f2py intent(in) nifx, nify, nf
    !f2py intent(in) kernel
    !f2py intent(in) iin

    ! input parameters
    integer :: ninx, niny
    integer :: nifx, nify, nf
    real*8, dimension(:, :) :: kernel
    real, dimension(:, :) :: iin

    ! auxiliar variables
    integer, target :: k, l

    iout(1: nifx, 1: nify) = 0

    write (*,*) shape(Iin), shape(iout)

    ! filter loop
    !$OMP PARALLEL DO private (k, l) shared (nf, nifx, nify, kernel, Iin, iout)
    do l = 1, nf
       do k = 1, nf
          iout(1: nifx, 1: nify) &
               = iout(1: nifx, 1: nify) + kernel(k, l) * Iin(k: Nifx + k, l: Nify + l)
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine conv

!  subroutine conv2 (ninx, niny, Nifx, Nify, nf, kernel1, kernel2, i1in, i2in)
!
!    implicit none
!
!    ! input parameters
!    integer :: ninx, niny
!    integer :: nifx, nify, nf
!    real*8, dimension(:, :) :: kernel1, kernel2
!    real, dimension(:, :) :: i1in, i2in
!
!    ! auxiliar variables
!    integer, target :: k, l
!
!    i1out(1: nifx, 1: nify) = 0
!    i2out(1: nifx, 1: nify) = 0
!
!    ! filter loop
!    do k = 1, nf
!       do l = 1, nf
!          i1out(1: nifx, 1: nify) = i1out(1: nifx, 1: nify) + kernel1(k, l) * I1in(k: Nifx + k, l: Nify + l)
!          i2out(1: nifx, 1: nify) = i2out(1: nifx, 1: nify) + kernel2(k, l) * I2in(k: Nifx + k, l: Nify + l)
!       end do
!    end do
!
!  end subroutine conv2

end module conv2fast
