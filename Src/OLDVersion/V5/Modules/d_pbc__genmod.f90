        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 18 15:09:44 2013
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
