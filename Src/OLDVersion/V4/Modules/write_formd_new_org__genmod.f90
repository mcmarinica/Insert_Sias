        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 16 12:09:22 2010
        MODULE WRITE_FORMD_NEW_ORG__genmod
          INTERFACE 
            SUBROUTINE WRITE_FORMD_NEW_ORG(NTEMP,NSIZE,RX,RY,RZ,A,A0,   &
     &VOLC)
              INTEGER(KIND=4), INTENT(IN) :: NSIZE
              INTEGER(KIND=4), INTENT(IN) :: NTEMP
              REAL(KIND=8), INTENT(IN) :: RX(NSIZE)
              REAL(KIND=8), INTENT(IN) :: RY(NSIZE)
              REAL(KIND=8), INTENT(IN) :: RZ(NSIZE)
              REAL(KIND=8), INTENT(IN) :: A(3,3)
              REAL(KIND=8), INTENT(IN) :: A0
              REAL(KIND=8), INTENT(IN) :: VOLC
            END SUBROUTINE WRITE_FORMD_NEW_ORG
          END INTERFACE 
        END MODULE WRITE_FORMD_NEW_ORG__genmod
