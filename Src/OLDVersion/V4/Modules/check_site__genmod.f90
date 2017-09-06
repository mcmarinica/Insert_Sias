        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 16 12:09:23 2010
        MODULE CHECK_SITE__genmod
          INTERFACE 
            SUBROUTINE CHECK_SITE(TSIA,NTEMP,NSIZE,NSIA,N1,N2,N3,NS1,NS2&
     &,NS3,NSITE,NSITE_M,NSITE_OCTA,TYPE_OCTA)
              INTEGER(KIND=4) :: NSIA
              INTEGER(KIND=4) :: NSIZE
              INTEGER(KIND=4) :: TSIA
              INTEGER(KIND=4) :: NTEMP
              INTEGER(KIND=4) :: N1(NSIZE)
              INTEGER(KIND=4) :: N2(NSIZE)
              INTEGER(KIND=4) :: N3(NSIZE)
              INTEGER(KIND=4) :: NS1(NSIA)
              INTEGER(KIND=4) :: NS2(NSIA)
              INTEGER(KIND=4) :: NS3(NSIA)
              INTEGER(KIND=4) :: NSITE(NSIA)
              INTEGER(KIND=4) :: NSITE_M(4,NSIA)
              INTEGER(KIND=4) :: NSITE_OCTA(10,NSIA)
              INTEGER(KIND=4) :: TYPE_OCTA(NSIA)
            END SUBROUTINE CHECK_SITE
          END INTERFACE 
        END MODULE CHECK_SITE__genmod
