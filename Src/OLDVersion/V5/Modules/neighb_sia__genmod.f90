        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 18 15:09:44 2013
        MODULE NEIGHB_SIA__genmod
          INTERFACE 
            SUBROUTINE NEIGHB_SIA(NALL,NSIA,VOLC,NSITE,RX,RY,RZ,A)
              INTEGER(KIND=4) :: NSIA
              INTEGER(KIND=4) :: NALL
              REAL(KIND=8) :: VOLC
              INTEGER(KIND=4) :: NSITE(NSIA)
              REAL(KIND=8) :: RX(NALL)
              REAL(KIND=8) :: RY(NALL)
              REAL(KIND=8) :: RZ(NALL)
              REAL(KIND=8) :: A(3,3)
            END SUBROUTINE NEIGHB_SIA
          END INTERFACE 
        END MODULE NEIGHB_SIA__genmod
