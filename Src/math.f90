MODULE math

CONTAINS

  !-----------------------------

  FUNCTION matdet(A) result(det)

    implicit none

    REAL(kind(0.d0)), dimension(:,:), intent(in) :: A
    REAL(kind(0.d0)) :: det

    IF ( (size(A,1).NE.3).AND.(size(A,2).NE.3) ) &
        STOP '< MatDet >: size of matrix to invert should be equal to 3'
    
    det =  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
         + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1) &
         - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3)

  END FUNCTION matdet

  !-------------------------------

  SUBROUTINE matinv(A, B)

    implicit none

    REAL(kind(0.d0)), dimension(:,:), intent(in) :: A
    REAL(kind(0.d0)), dimension(:,:), intent(out) :: B

    REAL(kind(0.d0)) :: invdet

    ! Variables used with Lapack subroutine
    INTEGER :: i, n, INFO
    INTEGER, dimension(:), allocatable :: IPIV
    REAL(kind(0.d0)), dimension(:,:), allocatable :: tempA, invA

    if (size(A,1).NE.size(A,2)) &
         STOP '< MatInv >: matrix to invert is not square'

    SELECT CASE(size(A,1))
    CASE(1)
       B(1,1)=1.d0/A(1,1)
       
    CASE(2)
       invdet=1.d0/( A(1,1)*A(2,2)-A(1,2)*A(2,1) )
       B(1,1)=A(2,2)*invdet
       B(2,2)=A(1,1)*invdet
       B(1,2)=-A(2,1)*invdet
       B(2,1)=-A(1,2)*invdet
       
    CASE(3)
       invdet=1.d0/matdet(A)
       
       b(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
       b(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
       b(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
       
       b(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
       b(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
       b(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)
       
       b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
       b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
       b(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

       b(1:3,1:3)=b(1:3,1:3)*invdet
       
    CASE DEFAULT
!!$         STOP '< MatInv >: size of matrix to invert should be less than 3 &
!!$            &(lapack not implemented)'
       ! Use Lapack subroutine DGESV       
       n=size(A,1)
       allocate(tempA(n,n))
       tempA(1:n,1:n)=A(1:n,1:n) ! Needed in order to not modify A
       allocate(IPIV(1:n))
       allocate(invA(n,n))
       invA(1:n,1:n)=0.d0
       do i=1,n
          invA(i,i)=1.d0
       enddo
       call DGESV(n, n, tempA, n, IPIV, invA, n, INFO)
       if (INFO.NE.0) then
          write(0,*) 'Result of DGESV subroutine: INFO=', INFO 
          write(0,*)' < 0: if INFO = -i, the i-th argument had an illegal value'
          write(0,*)' > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization'
          write(0,*)'has been completed, but the factor U is exactly'
          write(0,*)'singular, so the solution could not be computed.'
          STOP'< MatInv >'
       endif
       B(1:n,1:n)=invA(1:n,1:n)
       deallocate(tempA, invA, IPIV)

    END SELECT

  END SUBROUTINE matinv
  
  !-------------------------------------------------------------------
  FUNCTION VectProd(a1,a2) RESULT(b)

    IMPLICIT NONE

    REAL(kind(0.d0)), dimension(1:3), intent(in) :: a1, a2
    REAL(kind(0.d0)), dimension(1:3) :: b

    b(1)=a1(2)*a2(3)-a1(3)*a2(2)
    b(2)=a1(3)*a2(1)-a1(1)*a2(3)
    b(3)=a1(1)*a2(2)-a1(2)*a2(1)

  END FUNCTION VectProd

 subroutine gcd(na,nb,ngcd)
    implicit none
    integer, intent (in) :: na,nb
    integer, intent (out) :: ngcd
    integer :: ia,ib,itemp

    ia=na
    ib=nb
    
    do while  (ib /= 0)
     itemp = ia
     ia = ib
     ib =mod(itemp,ib)
    end do
    ngcd=ia
   return
 end subroutine gcd







END MODULE math
