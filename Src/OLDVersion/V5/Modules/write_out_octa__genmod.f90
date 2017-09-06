        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 18 15:12:12 2013
        MODULE WRITE_OUT_OCTA__genmod
          INTERFACE 
            SUBROUTINE WRITE_OUT_OCTA(NALL,NTEMP,NEW_SOLID,GNDATOMS,VOLC&
     &,A0,RX,RY,RZ,RXB,RYB,RZB,A,NDATOMS,ATOMS_TO_REMOVE,               &
     &ATOMS_TO_REMOVE_TYPE)
              INTEGER(KIND=4), INTENT(IN) :: NDATOMS
              INTEGER(KIND=4) :: NALL
              INTEGER(KIND=4) :: NTEMP
              INTEGER(KIND=4) :: NEW_SOLID
              INTEGER(KIND=4) :: GNDATOMS
              REAL(KIND=8) :: VOLC
              REAL(KIND=8) :: A0
              REAL(KIND=8) :: RX(NALL)
              REAL(KIND=8) :: RY(NALL)
              REAL(KIND=8) :: RZ(NALL)
              REAL(KIND=8), INTENT(IN) :: RXB(NALL)
              REAL(KIND=8), INTENT(IN) :: RYB(NALL)
              REAL(KIND=8), INTENT(IN) :: RZB(NALL)
              REAL(KIND=8) :: A(3,3)
              REAL(KIND=8), INTENT(IN) :: ATOMS_TO_REMOVE(NDATOMS)
              INTEGER(KIND=4), INTENT(IN) :: ATOMS_TO_REMOVE_TYPE(      &
     &NDATOMS)
            END SUBROUTINE WRITE_OUT_OCTA
          END INTERFACE 
        END MODULE WRITE_OUT_OCTA__genmod
