MODULE common_data

  SAVE

  ! Lattice vector coordinates (A): at(1:3,1), at(2:3,2), ...
  ! and inverse matrix
  REAL(kind(0.d0)), dimension(1:3, 1:3) :: at, inv_at
  ! Maximal number of atoms in simulation box (used to allocate tables)
  INTEGER :: imm=0
  ! Real number of atoms in simulation box
  INTEGER :: im=0,im0=0
  ! Atom reduced coordinates in simulation box: xc(1:3,:)
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xc,xc0
   ! Atom real coordinates: xp(1:3,:)
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp,xp0
  ! Atom type
  INTEGER, dimension(:), allocatable :: ityp,ityp0
  ! the structure 1-> FCC ; 2-> BCC
  INTEGER                  :: bulk_structure=2
  !1 for RadialAnalisys, 0 put the R1 by hand.
  INTEGER                  :: RadialAnalysis=1
  INTEGER                  :: max_nVoisins_WS
  REAL(kind(0.d0)) :: Rinst
  REAL(kind(0.d0)), parameter :: low_limit=10.d0*epsilon(1.d0)
  REAL(kind(0.d0)):: zero=0.d0
END MODULE common_data
