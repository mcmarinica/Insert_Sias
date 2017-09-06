        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 18 15:12:12 2013
        MODULE WRITE_COORDINATES__genmod
          INTERFACE 
            SUBROUTINE WRITE_COORDINATES(NTEMP,NSIZE,NEW_SOLID,RXGAOS,  &
     &RYGAOS,RZGAOS,GNDATOMS,ATOMS_TO_REMOVE,NDATOMS,RX,RY,RZ)
              INTEGER(KIND=4), INTENT(IN) :: NDATOMS
              INTEGER(KIND=4), INTENT(IN) :: GNDATOMS
              INTEGER(KIND=4), INTENT(IN) :: NSIZE
              INTEGER(KIND=4), INTENT(IN) :: NTEMP
              INTEGER(KIND=4), INTENT(OUT) :: NEW_SOLID
              REAL(KIND=8), INTENT(IN) :: RXGAOS(GNDATOMS)
              REAL(KIND=8), INTENT(IN) :: RYGAOS(GNDATOMS)
              REAL(KIND=8), INTENT(IN) :: RZGAOS(GNDATOMS)
              REAL(KIND=8), INTENT(IN) :: ATOMS_TO_REMOVE(NDATOMS)
              REAL(KIND=8), INTENT(INOUT) :: RX(NSIZE)
              REAL(KIND=8), INTENT(INOUT) :: RY(NSIZE)
              REAL(KIND=8), INTENT(INOUT) :: RZ(NSIZE)
            END SUBROUTINE WRITE_COORDINATES
          END INTERFACE 
        END MODULE WRITE_COORDINATES__genmod
