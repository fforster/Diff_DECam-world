! fast convolution routine using fortran 90

module projfast
  
  implicit none
  include "omp_lib.h"
  
  real, dimension(4100, 4100) :: imageout
  real, parameter :: pi = 3.14159265359
  real, parameter :: invpi = 0.31830988618
  
contains
  
  ! set number of threads
  subroutine set_num_threads(n)

    implicit none

    !f2py threadsafe
    !f2py intent(in) n

    INTEGER :: n
    CALL OMP_SET_NUM_THREADS(n)
    
  end subroutine set_num_threads

  ! sinc function (digital processing version)
  real function sinc (x)
    
    implicit none

    real :: x

    if (abs(x) < 1.e-8) then
       sinc = 1.
    else
       sinc = invpi * sin(pi * x) / x 
    endif
    
  end function sinc

  ! linear transformation + Lanczos interpolation
  subroutine o0_lanczos (alanczos, nx1, ny1, nx2, ny2, a11, a22, b1, b2, imagein)
    
    implicit none
    
    !f2py threadsafe
    !f2py intent(in) alanczos
    !f2py intent(in) nx1, ny1, nx2, ny2
    !f2py intent(in) a11, a22, b1, b2
    !f2py intent(in) imagein

    ! input parameters
    integer, intent(in) :: alanczos
    integer, intent(in) :: nx1, ny1, nx2, ny2
    real, intent(in) :: a11, a22, b1, b2 
    real, intent(in), dimension(:, :) :: imagein

    ! auxiliar variables
    integer :: i, j, k, l
    integer :: xf, yf
    integer :: xc, yc
    real :: xi, yi
    real :: xl, Lxl, yl, Lyl, Lxl2Lyl2

    write (*, *) nx1, ny1, nx2, ny2

    ! initialize output image
    imageout = 0

    ! Lanczos interpolation
    !$OMP PARALLEL DO private (k, l, i, j, xi, yi, xf, xc, yf, yc, xl, Lxl, yl, Lyl, Lxl2Lyl2) shared(a11, a22, b1, b2, alanczos, nx1, ny1, nx2, ny2, imagein, imageout)
    do k = 1 - alanczos, alanczos
       do l = 1 - alanczos, alanczos
          do j = 1, ny1
             do i = 1, nx1

                ! transformed pixels
                xi = a11 * i + b1
                yi = a22 * j + b2
       
                ! we are inside the definition range
                if (xi >= alanczos .and. xi <= nx2 - alanczos .and. yi >= alanczos .and. yi <= ny2 - alanczos) then
                   
                   xf = int(floor(xi))
                   xc = int(ceiling(xi))
                   yf = int(floor(yi))
                   yc = int(ceiling(yi))
                   
                   xl = xi - (xf + k)
                   Lxl = sinc(xl) * sinc(xl / alanczos)
                   yl = yi - (yf + l)
                   Lyl = sinc(yl) * sinc(yl / alanczos)

                   imageout(i, j) = imageout(i, j) + imagein(xf + k, yf + l) * Lxl * Lyl

                end if
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    write (*,*) maxval(imageout)

  end subroutine o0_lanczos

  ! linear transformation + Lanczos interpolation
  subroutine o1_lanczos (alanczos, nx1, ny1, nx2, ny2, a11, a12, a21, a22, b1, b2, imagein)
    
    implicit none
    
    !f2py threadsafe
    !f2py intent(in) alanczos
    !f2py intent(in) nx1, ny1, nx2, ny2
    !f2py intent(in) a11, a12, a21, a22, b1, b2
    !f2py intent(in) imagein

    ! input parameters
    integer, intent(in) :: alanczos
    integer, intent(in) :: nx1, ny1, nx2, ny2
    real, intent(in) :: a11, a12, a21, a22, b1, b2 
    real, intent(in), dimension(:, :) :: imagein

    ! auxiliar variables
    integer :: i, j, k, l
    integer :: xf, yf
    integer :: xc, yc
    real :: xi, yi
    real :: xl, Lxl, yl, Lyl, Lxl2Lyl2

    write (*, *) nx1, ny1, nx2, ny2

    ! initialize output image
    imageout = 0

    ! Lanczos interpolation
    !$OMP PARALLEL DO private (k, l, i, j, xi, yi, xf, xc, yf, yc, xl, Lxl, yl, Lyl, Lxl2Lyl2) shared(a11, a12, a21, a22, b1, b2, alanczos, nx1, ny1, nx2, ny2, imagein, imageout)
    do k = 1 - alanczos, alanczos
       do l = 1 - alanczos, alanczos
          do j = 1, ny1
             do i = 1, nx1

                ! transformed pixels
                xi = a11 * i + a12 * j + b1
                yi = a21 * i + a22 * j + b2
       
                ! we are inside the definition range
                if (xi >= alanczos .and. xi <= nx2 - alanczos .and. yi >= alanczos .and. yi <= ny2 - alanczos) then
                   
                   xf = int(floor(xi))
                   xc = int(ceiling(xi))
                   yf = int(floor(yi))
                   yc = int(ceiling(yi))
                   
                   xl = xi - (xf + k)
                   Lxl = sinc(xl) * sinc(xl / alanczos)
                   yl = yi - (yf + l)
                   Lyl = sinc(yl) * sinc(yl / alanczos)

                   imageout(i, j) = imageout(i, j) + imagein(xf + k, yf + l) * Lxl * Lyl

                end if
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    write (*,*) maxval(imageout)

  end subroutine o1_lanczos

  ! order 2 transformation + Lanczos interpolation
  subroutine o2_lanczos (alanczos, nx1, ny1, nx2, ny2, a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22, imagein)
    
    implicit none
    
    !f2py threadsafe
    !f2py intent(in) alanczos
    !f2py intent(in) nx1, ny1, nx2, ny2
    !f2py intent(in) a11, a12, a21, a22, b1, b2
    !f2py intent(in) imagein

    ! input parameters
    integer, intent(in) :: alanczos
    integer, intent(in) :: nx1, ny1, nx2, ny2
    real, intent(in) :: a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22
    real, intent(in), dimension(:, :) :: imagein

    ! auxiliar variables
    integer :: i, j, k, l
    integer :: xf, yf
    integer :: xc, yc
    real :: xi, yi
    real :: xl, Lxl, yl, Lyl, Lxl2Lyl2

    write (*, *) nx1, ny1, nx2, ny2

    ! initialize output image
    imageout = 0

    ! Lanczos interpolation
    !$OMP PARALLEL DO private (k, l, i, j, xi, yi, xf, xc, yf, yc, xl, Lxl, yl, Lyl, Lxl2Lyl2) shared(a1, b11, b12, c1, d11, d12, a2, b21, b22, c2, d21, d22, alanczos, nx1, ny1, nx2, ny2, imagein, imageout)
    do k = 1 - alanczos, alanczos
       do l = 1 - alanczos, alanczos
          do j = 1, ny1
             do i = 1, nx1

                ! transformed pixels
                xi = a1 + b11 * i + b12 * j + c1 * i * j + d11 * i * i + d12 * j * j
                yi = a2 + b21 * i + b22 * j + c2 * i * j + d21 * i * i + d22 * j * j
       
                ! we are inside the definition range
                if (xi >= alanczos .and. xi <= nx2 - alanczos .and. yi >= alanczos .and. yi <= ny2 - alanczos) then
                   
                   xf = int(floor(xi))
                   xc = int(ceiling(xi))
                   yf = int(floor(yi))
                   yc = int(ceiling(yi))
                   
                   xl = xi - (xf + k)
                   Lxl = sinc(xl) * sinc(xl / alanczos)
                   yl = yi - (yf + l)
                   Lyl = sinc(yl) * sinc(yl / alanczos)

                   imageout(i, j) = imageout(i, j) + imagein(xf + k, yf + l) * Lxl * Lyl

                end if
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    write (*,*) maxval(imageout)

  end subroutine o2_lanczos

  ! order 3 transformation + Lanczos interpolation
  subroutine o3_lanczos (alanczos, nx1, ny1, nx2, ny2, a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22, imagein)
    
    implicit none
    
    !f2py threadsafe
    !f2py intent(in) alanczos
    !f2py intent(in) nx1, ny1, nx2, ny2
    !f2py intent(in) a11, a12, a21, a22, b1, b2
    !f2py intent(in) imagein

    ! input parameters
    integer, intent(in) :: alanczos
    integer, intent(in) :: nx1, ny1, nx2, ny2
    real, intent(in) :: a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22
    real, intent(in), dimension(:, :) :: imagein

    ! auxiliar variables
    integer :: i, j, k, l
    integer :: xf, yf
    integer :: xc, yc
    real :: xi, yi
    real :: xl, Lxl, yl, Lyl, Lxl2Lyl2

    write (*, *) nx1, ny1, nx2, ny2

    ! initialize output image
    imageout = 0

    ! Lanczos interpolation
    !$OMP PARALLEL DO private (k, l, i, j, xi, yi, xf, xc, yf, yc, xl, Lxl, yl, Lyl, Lxl2Lyl2) shared(a1, b11, b12, c1, d11, d12, e1, f1, g11, g12, a2, b21, b22, c2, d21, d22, e2, f2, g21, g22, alanczos, nx1, ny1, nx2, ny2, imagein, imageout)
    do k = 1 - alanczos, alanczos
       do l = 1 - alanczos, alanczos
          do j = 1, ny1
             do i = 1, nx1

                ! transformed pixels
                xi = a1 + b11 * i + b12 * j + c1 * i * j + d11 * i * i + d12 * j * j + e1 * i * i * j + f1 * i * j * j + g11 * i * i * i + g12 * j * j * j
                yi = a2 + b21 * i + b22 * j + c2 * i * j + d21 * i * i + d22 * j * j + e2 * i * i * j + f2 * i * j * j + g21 * i * i * i + g22 * j * j * j
       
                ! we are inside the definition range
                if (xi >= alanczos .and. xi <= nx2 - alanczos .and. yi >= alanczos .and. yi <= ny2 - alanczos) then
                   
                   xf = int(floor(xi))
                   xc = int(ceiling(xi))
                   yf = int(floor(yi))
                   yc = int(ceiling(yi))
                   
                   xl = xi - (xf + k)
                   Lxl = sinc(xl) * sinc(xl / alanczos)
                   yl = yi - (yf + l)
                   Lyl = sinc(yl) * sinc(yl / alanczos)

                   imageout(i, j) = imageout(i, j) + imagein(xf + k, yf + l) * Lxl * Lyl

                end if
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    write (*,*) maxval(imageout)

  end subroutine o3_lanczos

end module projfast
