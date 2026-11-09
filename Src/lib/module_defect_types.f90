!=======================================================================
! module_defect_types.f90
!
! Module defect_types: the defect types of put_vct.x, as named constants.
! The value of each constant is the number typed at the prompt of the
! program, so a test reads  tsia==DEFECT_VACANCY  instead of  tsia==5.
!
! The functions group the types that the program treats alike:
!
!   defect_is_known      any type of the list
!   defect_is_dumbbell   dumbbells put on the listed sites
!   defect_is_migration  dumbbell jumps
!   defect_is_c15        C15 clusters
!   defect_is_ico        icosahedral clusters
!   defect_is_cluster    C15 or icosahedral clusters
!=======================================================================
module defect_types

  implicit none
  public

  ! dumbbells: one on each listed site
  integer, parameter :: DEFECT_SIA_110      = 1     ! <110> dumbbell
  integer, parameter :: DEFECT_SIA_111      = 2     ! <111> dumbbell
  integer, parameter :: DEFECT_SIA_100      = 3     ! <100> dumbbell
  integer, parameter :: DEFECT_SIA_111_CARE = 21    ! <111> dumbbell, neighbours of the row displaced
  integer, parameter :: DEFECT_SIA_100_CARE = 31    ! <100> dumbbell, neighbours of the row displaced
  integer, parameter :: DEFECT_FRANK_LOOP   = 210   ! <111> dumbbell with one atom shifted, for Frank loops in fcc

  ! dumbbell jumps: initial and final states of a <111> dumbbell
  integer, parameter :: DEFECT_SIA_111_MIGRATION      = 22
  integer, parameter :: DEFECT_SIA_111_MIGRATION_CARE = 221   ! with the neighbours of the row displaced

  ! other point defects
  integer, parameter :: DEFECT_GAO     = 4          ! Gao di-interstitial, from gao_unit.in
  integer, parameter :: DEFECT_VACANCY = 5          ! vacancy

  ! clusters
  integer, parameter :: DEFECT_C15           = 6    ! C15 cluster
  integer, parameter :: DEFECT_C15_TWISTED   = 67   ! C15 cluster, also cut out in a rotated cell
  integer, parameter :: DEFECT_ICO           = 7    ! icosahedral cluster
  integer, parameter :: DEFECT_ICO_DUMBBELLS = 77   ! icosahedral cluster with extra dumbbells and atoms

contains

!-----------------------------------------------------------------------
  pure logical function defect_is_known (tsia)
    integer, intent(in) :: tsia

    defect_is_known = defect_is_dumbbell(tsia) .or. defect_is_migration(tsia) .or. &
                      defect_is_cluster(tsia)  .or. (tsia==DEFECT_GAO)        .or. &
                      (tsia==DEFECT_VACANCY)

  end function defect_is_known

!-----------------------------------------------------------------------
  pure logical function defect_is_dumbbell (tsia)
    integer, intent(in) :: tsia

    defect_is_dumbbell = (tsia==DEFECT_SIA_110)      .or. (tsia==DEFECT_SIA_111)      .or. &
                         (tsia==DEFECT_SIA_100)      .or. (tsia==DEFECT_SIA_111_CARE) .or. &
                         (tsia==DEFECT_SIA_100_CARE) .or. (tsia==DEFECT_FRANK_LOOP)

  end function defect_is_dumbbell

!-----------------------------------------------------------------------
  pure logical function defect_is_migration (tsia)
    integer, intent(in) :: tsia

    defect_is_migration = (tsia==DEFECT_SIA_111_MIGRATION) .or. &
                          (tsia==DEFECT_SIA_111_MIGRATION_CARE)

  end function defect_is_migration

!-----------------------------------------------------------------------
  pure logical function defect_is_c15 (tsia)
    integer, intent(in) :: tsia

    defect_is_c15 = (tsia==DEFECT_C15) .or. (tsia==DEFECT_C15_TWISTED)

  end function defect_is_c15

!-----------------------------------------------------------------------
  pure logical function defect_is_ico (tsia)
    integer, intent(in) :: tsia

    defect_is_ico = (tsia==DEFECT_ICO) .or. (tsia==DEFECT_ICO_DUMBBELLS)

  end function defect_is_ico

!-----------------------------------------------------------------------
  pure logical function defect_is_cluster (tsia)
    integer, intent(in) :: tsia

    defect_is_cluster = defect_is_c15(tsia) .or. defect_is_ico(tsia)

  end function defect_is_cluster

end module defect_types
