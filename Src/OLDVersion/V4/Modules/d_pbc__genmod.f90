        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 16 12:09:23 2010
        MODULE D_PBC__genmod
          INTERFACE 
            SUBROUTINE D_PBC(UNIT_VECT,TAU,DS,NAT_UP)
              INTEGER(KIND=4) :: NAT_UP
              REAL(KIND=8) :: UNIT_VECT(3,3)
              REAL(KIND=8) :: TAU(3,NAT_UP)
              REAL(KIND=8) :: DS(NAT_UP*NAT_UP)
            END SUBROUTINE D_PBC
          END INTERFACE 
        END MODULE D_PBC__genmod
