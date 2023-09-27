    SUBROUTINE H2O2POT(v, x)
    IMPLICIT NONE
      REAL*8, INTENT(OUT) :: v
      REAL*8, INTENT(IN) :: x(12)
      real*8 :: q1, q2, q3, q4, q5, q6
      real*8 :: r, theta
      real*8,parameter :: Pi=acos(-1d0)

      r = sqrt(sum((x(4:6)-x(7:9))**2))
      q1 = 1d0 - 1.45625d0/r
      r = sqrt(sum((x(4:6)-x(1:3))**2))
      q2 = 1d0 - 0.96272d0/r
      r = sqrt(sum((x(7:9)-x(10:12))**2))
      q3 = 1d0 - 0.96272d0/r
      call thetaijk(x(1:3), x(4:6), x(6:9), theta)
      q4 = theta-101.034d0
      call thetaijk(x(4:6), x(6:9), x(10:12), theta)
      q5 = theta-101.034d0
      call dihedral(x(1:3), x(4:6), x(6:9), x(10:12), theta)
      q6 = theta
      v = -151.40422461    &
      &+ 0.03057888 * q2**1*q4**3 &
      &+ 0.00482409 * cos(q6) &
      &+ -0.06363999 * q2**3*q4**1 &
      &+ 0.00325127 * cos(2 * q6) &
      &+ -0.00373964 * q4**1*q5**3 &
      &+ 0.00026630 * cos(3 * q6) &
      &+ -0.04114668 * q1**2*q2**1*q3**1 &
      &+ 0.00004135 * cos(4 * q6) &
      &+ 0.11249614 * q1**2*q2**1*q4**1 &
      &+ 1.08154998 * q1**2 &
      &+ 0.02616679 * q1**2*q2**1*q5**1 &
      &+ 0.86097723 * q2**2 &
      &+ -0.07824425 * q1**2*q4**1*q5**1 &
      &+ 0.11016112 * q4**2 &
      &+ 0.04266205 * q1**1*q2**1*q4**2 &
      &+ -0.03637740 * q1**1*q2**1 &
      &+ -0.07420432 * q1**1*q2**2*q4**1 &
      &+ 0.17491854 * q1**1*q4**1 &
      &+ -0.08251268 * q1**1*q2**2*q5**1 &
      &+ 0.00057054 * q2**1*q3**1 &
      &+ 0.00270940 * q1**1*q4**1*q5**2 &
      &+ -0.00137967 * q2**1*q4**1 &
      &+ 0.00199953 * q2**1*q3**1*q4**2 &
      &+ 0.00062152 * q2**1*q5**1 &
      &+ -0.01292325 * q2**2*q4**1*q5**1 &
      &+ 0.01375726 * q4**1*q5**1 &
      &+ -0.02074323 * q2**1*q4**1*q5**2 &
      &+ -1.15691421 * q1**3 &
      &+ -0.00789732 * q2**1*q4**2*q5**1 &
      &+ -0.25495918 * q2**3 &
      &+ -0.01435326 * q1**1*cos(q6) &
      &+ -0.02272830 * q4**3 &
      &+ -0.00180710 * q2**1*cos(q6) &
      &+ -0.18381415 * q1**2*q2**1 &
      &+ -0.01135671 * q4**1*cos(q6) &
      &+ -0.34627237 * q1**2*q4**1 &
      &+ 0.00020655 * q1**2*cos(q6) &
      &+ 0.13974588 * q1**1*q2**2 &
      &+ -0.00492533 * q2**2*cos(q6) &
      &+ -0.27290455 * q1**1*q4**2 &
      &+ 0.00270990 * q4**2*cos(q6) &
      &+ -0.00674721 * q2**1*q3**2 &
      &+ 0.00376086 * q1**1*q2**1*cos(q6) &
      &+ -0.02179545 * q2**1*q4**2 &
      &+ 0.00044732 * q1**1*q4**1*cos(q6) &
      &+ -0.02125643 * q2**2*q4**1 &
      &+ 0.00569979 * q2**1*q3**1*cos(q6) &
      &+ -0.00491968 * q2**1*q5**2 &
      &+ -0.00244774 * q2**1*q4**1*cos(q6) &
      &+ 0.00233773 * q2**2*q5**1 &
      &+ -0.02065564 * q4**1*q5**1*cos(q6) &
      &+ -0.00050066 * q4**1*q5**2 &
      &+ 0.05249331 * q1**3*cos(q6) &
      &+ 0.01817536 * q1**1*q2**1*q3**1 &
      &+ -0.02490299 * q2**3*cos(q6) &
      &+ -0.04666009 * q1**1*q2**1*q4**1 &
      &+ 0.00391460 * q4**3*cos(q6) &
      &+ -0.02424748 * q1**1*q2**1*q5**1 &
      &+ 0.08893744 * q1**1*q2**2*cos(q6) &
      &+ -0.01727148 * q1**1*q4**1*q5**1 &
      &+ -0.01051618 * q1**1*cos(2 * q6) &
      &+ -0.00420506 * q2**1*q3**1*q4**1 &
      &+ 0.00120479 * q2**1*cos(2 * q6) &
      &+ -0.00647944 * q2**1*q4**1*q5**1 &
      &+ -0.00111888 * q4**1*cos(2 * q6) &
      &+ -1.06749007 * q1**4 &
      &+ 0.00884757 * q1**2*cos(2 * q6) &
      &+ -0.35741007 * q2**4 &
      &+ 0.00416289 * q2**2*cos(2 * q6) &
      &+ -0.00796836 * q4**4 &
      &+ 0.00126763 * q4**2*cos(2 * q6) &
      &+ -0.42556742 * q1**2*q2**2 &
      &+ -0.00706563 * q1**1*q2**1*cos(2 * q6) &
      &+ 0.06278896 * q1**2*q4**2 &
      &+ -0.00840146 * q1**1*q4**1*cos(2 * q6) &
      &+ -0.04010419 * q2**2*q4**2 &
      &+ -0.00139219 * q2**1*q3**1*cos(2 * q6) &
      &+ -0.00993912 * q4**2*q5**2 &
      &+ 0.00801673 * q2**1*q4**1*cos(2 * q6) &
      &+ 0.47562894 * q1**3*q2**1 &
      &+ 0.00463860 * q4**1*q5**1*cos(2 * q6) &
      &+ -0.40830627 * q1**3*q4**1 &
      &+ -0.00096051 * q1**1*cos(3 * q6) &
      &+ 0.22073222 * q1**1*q2**3 &
      &+ 0.00019906 * q2**1*cos(3 * q6) &
      &+ 0.07828212 * q1**1*q4**3 &
      &+ -0.00057576 * q4**1*cos(3 * q6) &
      &+ -0.02954687 * q2**1*q3**3
    END SUBROUTINE H2O2POT

    subroutine thetaijk(i, j, k, theta)
      implicit none
      real*8, intent(IN) :: i(3), j(3), k(3)
      real*8, intent(out) :: theta
      real*8 :: r1, r2, cos=0d0
      integer :: ii
      real*8,parameter :: Pi=acos(-1d0)  
      r1 = sqrt(sum((i-j)**2))
      r2 = sqrt(sum((k-j)**2))
      do ii = 1, 3
        cos = cos + (i(ii)-j(ii))*(k(ii)-j(ii))
      end do
      cos = cos/r1/r2
      theta = acos(cos) / Pi * 180d0
    end subroutine thetaijk
    
    function cross_product(vec1,vec2)
        real*8 :: cross_product(3)
        real*8 :: vec1(3),vec2(3)
        cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    end function cross_product
    
    subroutine crossproduct(vec1, vec2, prod)
        real*8,intent(in) :: vec1(3),vec2(3)
        real*8,intent(out) :: prod(3)
        prod(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        prod(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        prod(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    end subroutine crossproduct
    
    subroutine dihedral(a, b, c, d, theta)
    
        real*8,intent(in) :: a(3), b(3), c(3), d(3)
        real*8,intent(out) :: theta
        real*8 :: vec1(3), vec2(3)
        integer :: ii=0
        real*8, parameter :: Pi = acos(-1d0)
!        vec1(:) = cross_product(b-a, c-b)
!        vec2(:) = cross_product(c-b, d-c)
        call crossproduct(b-a, c-b, vec1)
        call crossproduct(c-b, d-c, vec2)
        theta = 0d0
        do ii=1, 3
            theta = theta + vec1(ii)*vec2(ii)
        end do
        theta = acos(theta/sqrt(sum(vec1**2)*sum(vec2**2)))
    end subroutine dihedral
        



