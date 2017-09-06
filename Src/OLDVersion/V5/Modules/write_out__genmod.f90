        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 18 15:09:44 2013
        MODULE WRITE_OUT__genmod
          INTERFACE 
            SUBROUTINE WRITE_OUT(NALL,NTEMP,NSIA,VOLC,NSITE,RX,RY,RZ,RSX&
     &,RSY,RSZ,A,A0)
              INTEGER(KIND=4) :: NSIA
              INTEGER(KIND=4) :: NALL
              INTEGER(KIND=4) :: NTEMP
              REAL(KIND=8) :: VOLC
              INTEGER(KIND=4) :: NSITE(NSIA)
              REAL(KIND=8) :: RX(NALL)
              REAL(KIND=8) :: RY(NALL)
              REAL(KIND=8) :: RZ(NALL)
              REAL(KIND=8) :: RSX(2)
              REAL(KIND=8) :: RSY(2)
              REAL(KIND=8) :: RSZ(2)
              REAL(KIND=8), INTENT(IN) :: A(3,3)
              REAL(KIND=8), INTENT(IN) :: A0
            END SUBROUTINE WRITE_OUT
          END INTERFACE 
        END MODULE WRITE_OUT__genmod
