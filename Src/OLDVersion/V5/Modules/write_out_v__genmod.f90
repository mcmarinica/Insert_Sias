        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 18 15:09:44 2013
        MODULE WRITE_OUT_V__genmod
          INTERFACE 
            SUBROUTINE WRITE_OUT_V(NALL,NTEMP,NSIA,NSIZE,VOLC,NSITE,RX, &
     &RY,RZ,A)
              INTEGER(KIND=4) :: NSIZE
              INTEGER(KIND=4) :: NSIA
              INTEGER(KIND=4) :: NALL
              INTEGER(KIND=4) :: NTEMP
              REAL(KIND=8) :: VOLC
              INTEGER(KIND=4) :: NSITE(NSIA)
              REAL(KIND=8) :: RX(NSIZE)
              REAL(KIND=8) :: RY(NSIZE)
              REAL(KIND=8) :: RZ(NSIZE)
              REAL(KIND=8) :: A(3,3)
            END SUBROUTINE WRITE_OUT_V
          END INTERFACE 
        END MODULE WRITE_OUT_V__genmod
