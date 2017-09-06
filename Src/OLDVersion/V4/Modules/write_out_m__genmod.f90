        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 16 12:09:23 2010
        MODULE WRITE_OUT_M__genmod
          INTERFACE 
            SUBROUTINE WRITE_OUT_M(NALL,NTEMP,NSIA,VOLC,A0,NSITE_M,RX,RY&
     &,RZ,RSX,RSY,RSZ,A)
              INTEGER(KIND=4) :: NSIA
              INTEGER(KIND=4) :: NALL
              INTEGER(KIND=4) :: NTEMP
              REAL(KIND=8) :: VOLC
              REAL(KIND=8) :: A0
              INTEGER(KIND=4) :: NSITE_M(4,NSIA)
              REAL(KIND=8) :: RX(NALL)
              REAL(KIND=8) :: RY(NALL)
              REAL(KIND=8) :: RZ(NALL)
              REAL(KIND=8) :: RSX(6)
              REAL(KIND=8) :: RSY(6)
              REAL(KIND=8) :: RSZ(6)
              REAL(KIND=8) :: A(3,3)
            END SUBROUTINE WRITE_OUT_M
          END INTERFACE 
        END MODULE WRITE_OUT_M__genmod
