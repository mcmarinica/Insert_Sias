!=======================================================================
! put_vct.f90
!
! Inserts defects in a perfect bcc or fcc crystal.
!
!   input    standard input  the defect type
!            formd.str       the reference crystal, written by cube.x
!            no00.inp        the lattice sites that receive a defect, in
!                            units of a0/2 (no00m.inp for dumbbell jumps)
!            gao_unit.in     the six atoms of a Gao di-interstitial
!
!   output   formd_new.str, formd_new_ndm_scl.str, formd_new.lmp, ...
!            the crystal with the defects; doc/files.md lists every file
!
! Steps:
!   1. read the defect type, the reference crystal and the list of sites
!   2. give every atom its integer coordinates, in units of a0/2
!   3. find the atoms that sit on the listed sites
!   4. build the defects and write the output files
!
! The defect type, tsia, is one of the constants of the module
! defect_types (lib/module_defect_types.f90): DEFECT_SIA_110,
! DEFECT_VACANCY, DEFECT_C15, ... The value of a constant is the number
! typed at the prompt.
!=======================================================================
program put_vct

  use math
  use octahedre
  use icosaedre
  use migration_sia
  use the_111_smooth
  use the_100_smooth
  use lammps_output
  use defect_types
  implicit none

  !---------------------------------------------------------------------
  ! Compiled parameters, in Angstrom.
  !   a0      lattice parameter of the reference crystal; formd.str must
  !           be built with this value
  !   as110   distance between the two atoms of a <110> dumbbell
  !   as111   the same for a <111> dumbbell
  !   as100   the same for a <100> dumbbell
  ! The first line holds the values for bcc. For loops in fcc, comment it
  ! out and activate the line marked fcc.
  !---------------------------------------------------------------------
  real(8), parameter :: a0=2.8553,as110=2.05d0,as111=1.90d0,as100=1.70d0
  !fcc real(8), parameter :: a0=2.8553,as110=1.55d0,as111=1.550d0,as100=1.550d0

  ! defect type and numbers of atoms
  integer :: tsia          ! defect type
  integer :: ntemp         ! atoms of the reference crystal
  integer :: nsia          ! entries of the list of sites
  integer :: nsia_a        ! atoms outside lattice sites, as estimated from nsia
  integer :: nall          ! upper bound of the final number of atoms
  integer :: nsize         ! dimension of the coordinate arrays
  integer :: new_solid     ! final number of atoms (cluster types)
  integer :: ndatoms       ! lattice atoms removed (cluster types)
  integer :: gndatoms      ! atoms inserted (cluster types)
  integer :: nlac_remove, natoms_extra, old_size
  integer :: n1MAX, n2MAX, n3MAX      ! largest integer coordinates of the crystal

  ! the crystal
  real(8), dimension(3,3) :: a                              ! cell vectors, one per row
  real(8) :: volc                                           ! scaling of the cell, kept at 1
  real(8), dimension(:), allocatable :: rx, ry, rz          ! positions, modified by the defects
  real(8), dimension(:), allocatable :: rxb, ryb, rzb       ! positions of the reference crystal
  integer, dimension(:), allocatable :: n1, n2, n3          ! integer coordinates of the atoms

  ! the list of sites
  integer, dimension(:), allocatable :: ns1, ns2, ns3       ! integer coordinates of the sites
  integer, dimension(:), allocatable :: nsite               ! atom that sits on each site
  integer, dimension(:,:), allocatable :: nsite_m           ! the four atoms a Gao di-interstitial removes
  integer, dimension(:), allocatable :: type_octa           ! type of each C15 unit
  integer, dimension(:), allocatable :: nd1, nd2, nd3       ! extra dumbbells and atoms of type 77

  ! dumbbells and Gao di-interstitials
  real(8) :: rsx(2), rsy(2), rsz(2)      ! the two atoms of a dumbbell, relative to the site
  real(8) :: rgx(6), rgy(6), rgz(6)      ! the six atoms of gao_unit.in

  ! clusters: removed lattice atoms and inserted atoms
  integer, dimension(:), allocatable :: nsite_octa_vac_remove, nsite_octa_vac_remove_type
  integer, dimension(:), allocatable :: indx, atoms_to_remove_type
  real(8), dimension(:), allocatable :: list_no_order, atoms_to_remove, offlattice_no_order
  real(8), dimension(:), allocatable :: rxgaos, rygaos, rzgaos
  real(8), dimension(:), allocatable :: temprxgaos, temprygaos, temprzgaos
  integer, dimension(:), allocatable :: buffer_atoms, buffer_type
  real(8), dimension(:), allocatable :: buffer_rx, buffer_ry, buffer_rz

  !=====================================================================
  ! 1. Input
  !=====================================================================

  call read_tsia (tsia)
  call read_ntemp (ntemp)          ! number of atoms of formd.str
  call read_nsia (nsia,tsia)       ! first line of no00.inp or no00m.inp

  ! atoms outside lattice sites and upper bound of the final number of atoms
  if (tsia==DEFECT_GAO) then
    nsia_a=nsia*4
    nall  =ntemp+2*nsia
  else if (tsia==DEFECT_VACANCY) then
    nsia_a=-nsia
    nall=ntemp-nsia
  else if (defect_is_c15(tsia)) then
    call init_octa
    nsia_a=nsia*10
    nall=ntemp + 2*nsia_a
  else if (tsia==DEFECT_ICO) then
    call init_icos
    nsia_a=13*nsia                           ! one atom in the centre and six dumbbells
    nall = ntemp+nsia_a
  else if (tsia==DEFECT_ICO_DUMBBELLS) then
    call init_icos
    nsia_a=13*nsia_ico + 2*nsia_dumbell
    nall = ntemp+nsia_a
  else
    nsia_a=nsia
    nall  =ntemp+nsia
  end if

  write(*,'("The number of lattice atoms .......:",i8)') ntemp
  write(*,'("The number of lines in no00.inp....:",i8)') nsia
  write(*,'("The number of no-lattice atoms ....:",i8)') nsia_a

  nsize=nall
  if (nall.lt.ntemp) nsize=ntemp
  allocate ( n1(nsize),  n2(nsize),  n3(nsize),                  &
             rx(nsize),  ry(nsize), rz(nsize),                   &
             rxb(nsize),  ryb(nsize), rzb(nsize))
  allocate (ns1(nsia),ns2(nsia),ns3(nsia),nsite(nsia),          &
            nsite_m(4,nsia),nsite_octa_vac_remove(10*nsia),     &
            nsite_octa_vac_remove_type(10*nsia),type_octa(nsia),&
            indx(10*nsia), list_no_order(10*nsia),              &
            offlattice_no_order(12*nsia), temprxgaos(12*nsia),  &
            temprygaos(12*nsia),temprzgaos(12*nsia))
  nsite(:)=0         ! 0: no atom sits on the site (centre of an icosahedral unit)
  if (tsia==DEFECT_ICO_DUMBBELLS) then
    allocate(nd1(nsia_dumbell), nd2(nsia_dumbell), nd3(nsia_dumbell), type_dumb(nsia_dumbell))
  end if

  ! the reference crystal: positions rx,ry,rz and cell a; the entries of
  ! the arrays beyond the ntemp atoms are ghosts, far from the crystal
  rx(:)=99999999.d0
  ry(:)=99999999.d0
  rz(:)=99999999.d0
  call read_input(ntemp,nsize,a,rx,ry,rz)
  call lammps_write_all ('formd.lmp', ntemp, a, rx, ry, rz, &
       'Insert_Sias: reference crystal (formd.str), lengths in Angstrom')

  ! geometry of the dumbbell or of the Gao di-interstitial
  if (nsia_a.gt.0) then
    call read_sia (tsia,as110,as111,as100,rsx,rsy,rsz,rgx,rgy,rgz)
    if (tsia==DEFECT_SIA_111_CARE) call allocate_the_111 (nsia)
    if (tsia==DEFECT_SIA_100_CARE) call allocate_the_100 (nsia)
  end if

  !=====================================================================
  ! 2. Integer coordinates of the atoms, in units of a0/2
  !=====================================================================

  call write_n123 (ntemp,nsize,a0,n1,n2,n3,rx,ry,rz,  &
                   n1MAX,n2MAX,n3MAX)

  !=====================================================================
  ! 3. The listed sites and the atoms that sit on them
  !=====================================================================

  ! every type except the dumbbell jumps: sites of no00.inp
  if (.not. defect_is_migration(tsia)) then

    call read_config_sia (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)

    ! orientation of the icosahedral units
    if (tsia==DEFECT_ICO)  call fix_ico_type (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
    if (tsia==DEFECT_ICO_DUMBBELLS) call fix_ico_type_mixed (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)

    ! neighbour statistics of the sites, fort.567
    if (tsia/=DEFECT_ICO_DUMBBELLS) then
      call analyze_config_sia (tsia,nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
    end if

    ! displacements of the row neighbours
    if (tsia==DEFECT_SIA_111_CARE) call fill_the_111 (a0,as111,rsx,rsy,rsz)
    if (tsia==DEFECT_SIA_100_CARE) call fill_the_100 (a0,as100,rsx,rsy,rsz)

    ! the rotated cut of type 67 needs a square cell of suitable size
    if (tsia==DEFECT_C15_TWISTED) then
      if ((n1MAX/=n2MAX).or.(mod(n1MAX-1,4)/=0)) then
        write(*,'("n1MAX-1 and n2MAX-1 should be equal and multiple of 4...:", 2i4)') n1MAX-1,n2MAX-1
        write(*,'("stop")')
        stop
      end if
    end if

    ! C15: type of each unit, one of the two orientations of the tetrahedron
    if (defect_is_c15(tsia)) then
      type_octa(:)=1
      call read_config_sia_octa (nsia,type_octa,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
    end if

  end if

  ! type 77: the extra dumbbells and single atoms
  if (tsia==DEFECT_ICO_DUMBBELLS) then
    call read_config_sia_ico (nsia_dumbell,nd1,nd2,nd3)
  end if

  ! dumbbell jumps: sites and jump flags of no00m.inp
  if (defect_is_migration(tsia)) then
    call allocate_migration(nsia)
    if (tsia==DEFECT_SIA_111_MIGRATION_CARE) then
      call allocate_the_111_migration (ntemp,nsia)
    end if
    call read_config_sia_migration (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
    if (tsia==DEFECT_SIA_111_MIGRATION_CARE)then
      call fill_the_111 (a0,as111,rsx,rsy,rsz)
    end if
  end if

  ! sites of the row neighbours
  if (tsia==DEFECT_SIA_111_CARE)  call build_the_nsite_the_111 (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
  if (tsia==DEFECT_SIA_100_CARE)  call build_the_nsite_the_100 (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
  if (tsia==DEFECT_SIA_111_MIGRATION_CARE) call build_the_nsite_the_111_migration (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)

  ! atom of each site, and the lattice atoms that the larger objects remove
  call check_site (tsia,ntemp,nsize,nsia,nlac_remove,n1,n2,n3,ns1,ns2,ns3,&
                   nsite,nsite_m,nsite_octa_vac_remove,nsite_octa_vac_remove_type,type_octa)
  if (defect_is_migration(tsia)) then
    call check_site_migration (ntemp,nsize,nsia,n1,n2,n3,ns1,ns2,ns3)
    if (tsia==DEFECT_SIA_111_MIGRATION_CARE) then
      call care_check_site_migration (ntemp,nsize,nsia,n1,n2,n3)
    end if
  end if

  !=====================================================================
  ! 4. The defects and the output files
  !=====================================================================

  volc=1.d0          ! the cell of the reference crystal is kept

  !---------------------------------------------------------------------
  ! dumbbells: types 1, 2, 3, 21, 31, 210
  !---------------------------------------------------------------------
  if (defect_is_dumbbell(tsia)) then

    rxb(:)=rx(:)
    ryb(:)=ry(:)
    rzb(:)=rz(:)

    if (tsia==DEFECT_SIA_111_CARE) then
      call distord_lattice_111 (ntemp,nsia,rx,ry,rz)
    end if
    if (tsia==DEFECT_SIA_100_CARE) then
      call distord_lattice_100 (ntemp,nsia,rx,ry,rz)
    end if
    call write_out   (nall,ntemp,nsia,volc,nsite,rxb,ryb,rzb,rx,ry, &
                      rz,rsx,rsy,rsz,a,a0)

  !---------------------------------------------------------------------
  ! Gao di-interstitials: type 4
  !---------------------------------------------------------------------
  else if (tsia==DEFECT_GAO) then

    call write_out_m (nall,ntemp,nsia,volc,a0,nsite_m,rx,ry, &
                      rz,rgx,rgy,rgz,a)

  !---------------------------------------------------------------------
  ! vacancies: type 5
  !---------------------------------------------------------------------
  else if (tsia==DEFECT_VACANCY) then

    call write_out_vv (nall,ntemp,nsia,nsize,volc,nsite,rx,ry,rz,a)

  !---------------------------------------------------------------------
  ! C15 and icosahedral clusters: types 6, 67, 7, 77
  !---------------------------------------------------------------------
  else if (defect_is_cluster(tsia)) then

    write(*,*) 'nlac_remove (before suppression).........:',nlac_remove

    call write_formd_new_org (ntemp,nsize,rx,ry,rz,a,a0,volc)
    rxb(:)=rx(:)
    ryb(:)=ry(:)
    rzb(:)=rz(:)

    ! Lattice atoms to remove. Linked units name some atoms twice:
    ! ndatoms is the number of distinct ones, atoms_to_remove their list.
    call  lattice_compute_different_removeatoms (nsia, nlac_remove,nsite_octa_vac_remove,indx,list_no_order,ndatoms)
    allocate (atoms_to_remove(ndatoms),atoms_to_remove_type(ndatoms))
    call  lattice_compute_atoms_to_remove(atoms_to_remove,atoms_to_remove_type,ndatoms, &
         nlac_remove,list_no_order,nsite_octa_vac_remove_type,indx)

    ! Atoms to insert, duplicates removed in the same way:
    ! gndatoms atoms at rxgaos,rygaos,rzgaos.
    if (defect_is_c15(tsia)) then

      ! twelve atoms per C15 unit
      if (allocated(indx)) deallocate(indx)             ; allocate(indx(12*nsia))

      call  offlattice_get_coordinates_and_order(a0,temprxgaos,temprygaos,temprzgaos,&
                                               rx,    ry,    rz,    &
                offlattice_no_order,indx,nsite,type_octa,nsia,nsize,gndatoms,natoms_extra)
      allocate (rxgaos(gndatoms),rygaos(gndatoms),rzgaos(gndatoms))
      call offlattice_compute_atoms_to_remove(gndatoms,rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos,&
                       offlattice_no_order,indx,nsia,natoms_extra)

    else if (defect_is_ico(tsia)) then

      ! thirteen atoms per icosahedral unit
      if (allocated(indx)) deallocate(indx)             ; allocate(indx(13*nsia))
      if (allocated(offlattice_no_order)) deallocate(offlattice_no_order) ; allocate(offlattice_no_order(13*nsia))
      if (allocated(temprxgaos)) deallocate(temprxgaos) ; allocate(temprxgaos(13*nsia))
      if (allocated(temprygaos)) deallocate(temprygaos) ; allocate(temprygaos(13*nsia))
      if (allocated(temprzgaos)) deallocate(temprzgaos) ; allocate(temprzgaos(13*nsia))

      call  ico_offlattice_get_coordinates_and_order(a0,temprxgaos,temprygaos,temprzgaos,&
                                   rx, ry, rz, ns1, ns2, ns3,   &
                offlattice_no_order,indx,nsite,type_octa,nsia,nsize,gndatoms,natoms_extra)
      allocate (rxgaos(gndatoms),rygaos(gndatoms),rzgaos(gndatoms))
      call ico_offlattice_compute_atoms_to_remove(gndatoms,rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos,&
                       offlattice_no_order,indx,nsia,natoms_extra)

    end if

    write(*,*)'ntemp - original bulk.....................:', ntemp
    write(*,*)'gndatoms - offlatice to be inserted.......:', gndatoms
    write(*,*)'ndatoms  - atoms to be removed ...........:', ndatoms
    write(*,*)'nsize - allocate dimension for rx.........:', nsize

    ! type 77: the extra dumbbells and single atoms join the two lists
    if ((tsia==DEFECT_ICO_DUMBBELLS).and.(nsia_dumbell > 0)) then

      call dumbell_ico(rx,ry,rz, ntemp, nd1,nd2,nd3,size(nd1), a0/2.d0)

      ! the lattice atoms that the extra dumbbells replace
      if (allocated(buffer_atoms)) deallocate(buffer_atoms) ; allocate(buffer_atoms(size(atoms_to_remove)))
      old_size = size(atoms_to_remove)
      if (old_size /= ndatoms) then
        write(*,*) 'problem in counting the middle of SIAs'
        stop
      end if
      buffer_atoms = atoms_to_remove
      if (allocated(atoms_to_remove)) deallocate(atoms_to_remove) ; allocate(atoms_to_remove(ndatoms+nsia_dumbell))
      atoms_to_remove(1:old_size) = buffer_atoms(1:old_size)
      atoms_to_remove(old_size+1:ndatoms+nsia_dumbell) = nsite_dumb_remove(1:nsia_dumbell)

      ! their type is 0: these sites get no marker in the visualisation files
      allocate(buffer_type(old_size))
      buffer_type = atoms_to_remove_type
      deallocate(atoms_to_remove_type) ; allocate(atoms_to_remove_type(ndatoms+nsia_dumbell))
      atoms_to_remove_type(1:old_size) = buffer_type(1:old_size)
      atoms_to_remove_type(old_size+1:ndatoms+nsia_dumbell) = 0

      ! the atoms of the extra dumbbells and the single atoms
      old_size = size(rxgaos)
      if (allocated(buffer_rx)) deallocate(buffer_rx) ; allocate(buffer_rx(old_size))
      if (allocated(buffer_ry)) deallocate(buffer_ry) ; allocate(buffer_ry(old_size))
      if (allocated(buffer_rz)) deallocate(buffer_rz) ; allocate(buffer_rz(old_size))
      buffer_rx = rxgaos
      buffer_ry = rygaos
      buffer_rz = rzgaos
      if (allocated(rxgaos)) deallocate(rxgaos) ; allocate(rxgaos(gndatoms + 2*nsia_dumbell_real_ico + nsia_center_real_ico))
      if (allocated(rygaos)) deallocate(rygaos) ; allocate(rygaos(gndatoms + 2*nsia_dumbell_real_ico + nsia_center_real_ico))
      if (allocated(rzgaos)) deallocate(rzgaos) ; allocate(rzgaos(gndatoms + 2*nsia_dumbell_real_ico + nsia_center_real_ico))

      rxgaos(1:gndatoms) = buffer_rx(1:gndatoms)
      rygaos(1:gndatoms) = buffer_ry(1:gndatoms)
      rzgaos(1:gndatoms) = buffer_rz(1:gndatoms)

      rxgaos(gndatoms+1:gndatoms+2*nsia_dumbell_real_ico + nsia_center_real_ico) = &
           r_extra_ico(1,1:2*nsia_dumbell_real_ico + nsia_center_real_ico)
      rygaos(gndatoms+1:gndatoms+2*nsia_dumbell_real_ico + nsia_center_real_ico) = &
           r_extra_ico(2,1:2*nsia_dumbell_real_ico + nsia_center_real_ico)
      rzgaos(gndatoms+1:gndatoms+2*nsia_dumbell_real_ico + nsia_center_real_ico) = &
           r_extra_ico(3,1:2*nsia_dumbell_real_ico + nsia_center_real_ico)

      gndatoms = gndatoms + 2*nsia_dumbell_real_ico + nsia_center_real_ico
      ndatoms  = ndatoms + nsia_dumbell

    end if

    ! final coordinates: the inserted atoms first, then the remaining crystal
    call write_coordinates (ntemp,nsize,new_solid,rxgaos,rygaos,rzgaos,gndatoms,atoms_to_remove,ndatoms,&
                            rx,ry,rz)

    write(*,*)'new_solid.................................:', new_solid
    write(*,*)'The no of defects SIAs ...................:', new_solid-ntemp

    call  write_out_octa (nall,ntemp,new_solid,gndatoms,volc,a0,rx,  &
                   ry,rz,rxb,ryb,rzb,a,ndatoms,atoms_to_remove,      &
                   atoms_to_remove_type, nsite, nsia)

    ! type 67: the cluster and the crystal cut out in a rotated cell
    if (tsia==DEFECT_C15_TWISTED) then
      call write_out_kazuto (a0,n1MAX,n3MAX,ntemp,new_solid,gndatoms,nsize,nsia,nsite, &
                             rx,ry,rz,rxb,ryb,rzb)
    end if

  !---------------------------------------------------------------------
  ! dumbbell jumps: types 22, 221
  !---------------------------------------------------------------------
  else if (defect_is_migration(tsia)) then

    rxb(:)=rx(:)
    ryb(:)=ry(:)
    rzb(:)=rz(:)

    ! initial state
    if (tsia==DEFECT_SIA_111_MIGRATION_CARE) then
      call deb_distord_lattice_111_migration (ntemp,nsia,rxb,ryb,rzb,rx,ry,rz)
    end if
    call deb_fill_the_mig_objects (ntemp,rx,ry,rz,rsx,rsy,rsz)
    call deb_write_out_migration   (nall,ntemp,volc,rx,ry, &
                    rz,rsx,rsy,rsz,a)

    call write_out   (nall,ntemp,nsia,volc,nsite,rxb,ryb,rzb,rx,ry, &
                    rz,rsx,rsy,rsz,a,a0)

    ! final state
    if (tsia==DEFECT_SIA_111_MIGRATION_CARE) then
      call fin_distord_lattice_111_migration (ntemp,nsia,rxb,ryb,rzb,rx,ry,rz)
    end if
    call fin_fill_the_mig_objects (ntemp,rx,ry,rz,rsx,rsy,rsz)
    call fin_write_out_migration   (nall,ntemp,volc,rx,ry, &
                    rz,rsx,rsy,rsz,a)

  end if

end program put_vct
