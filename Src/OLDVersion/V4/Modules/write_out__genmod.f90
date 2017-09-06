        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 16 12:09:23 2010
        MODULE WRITE_OUT__genmod
          INTERFACE 
            SUBROUTINE WRITE_OUT(NALL,NTEMP,NSIA,VOLC,NSITE,RX,RY,RZ,RSX&
     &,RSY,RSZ,A)
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
              REAL(KIND=8) :: A(3,3)
            END SUBROUTINE WRITE_OUT
          END INTERFACE 
        END MODULE WRITE_OUT__genmod
