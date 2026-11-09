!=======================================================================
! migration.f90
!
! Dumbbell jumps (types 22 and 221): the atoms of a jumping dumbbell in the
! initial and final states, and the files deb_loop.gin_0 and fin_loop.gin_0.
!
!   deb_fill_the_mig_objects  atoms of the jumping dumbbells, initial state
!   fin_fill_the_mig_objects  atoms of the jumping dumbbells, final state
!   deb_write_out_migration   writes the initial state
!   fin_write_out_migration   writes the final state
!=======================================================================

subroutine deb_fill_the_mig_objects (ntemp,rx,ry,rz,rsx,rsy,rsz)

use migration_sia
implicit none
integer :: ntemp
real(8) :: rx(ntemp),ry(ntemp),rz(ntemp)
real(8) :: rsx(2),rsy(2),rsz(2)

integer :: is,ii,jj


   do is=1,nmig
     ii = nsite_migra(is)
     jj = nsite_trans(is)
     debatoms(is)%atom(1:3,1) = (/ rx(ii)+rsx(1), ry(ii)+rsy(1), rz(ii)+rsz(1) /)   !first SIA
     debatoms(is)%atom(1:3,2) = (/ rx(ii)+rsx(2), ry(ii)+rsy(2), rz(ii)+rsz(2) /)   !second SIA
     debatoms(is)%atom(1:3,3) = (/ rx(jj), ry(jj),rz(jj) /) !<111> translation, bulk position of SIA
   end do
  return
end subroutine deb_fill_the_mig_objects

!=======================================================================

subroutine fin_fill_the_mig_objects (ntemp,rx,ry,rz,rsx,rsy,rsz)

use migration_sia
implicit none
integer :: ntemp
real(8) :: rx(ntemp),ry(ntemp),rz(ntemp)
real(8) :: rsx(2),rsy(2),rsz(2)

integer :: is,ii,jj


   do is=1,nmig
     ii = nsite_migra(is)
     jj = nsite_trans(is)
     finatoms(is)%atom(1:3,2) = (/ rx(jj)+rsx(1), ry(jj)+rsy(1), rz(jj)+rsz(1) /)   !first SIA
     finatoms(is)%atom(1:3,3) = (/ rx(jj)+rsx(2), ry(jj)+rsy(2), rz(jj)+rsz(2) /)   !second SIA
     finatoms(is)%atom(1:3,1) = (/ rx(ii), ry(ii),rz(ii) /) !<111> translation, bulk position of SIA
   end do
  return
end subroutine fin_fill_the_mig_objects

!=======================================================================

   subroutine deb_write_out_migration (nall,ntemp,volc,rx,  &
                  ry,rz,rsx,rsy,rsz,a)
  use math
  use pbc
  use migration_sia
  use the_111_smooth
  use lammps_output

  implicit none
  integer  :: nall,ntemp,i,j,ii,icard,is
  double precision,dimension(nall)     :: rx,ry,rz
  integer  :: in ,icnt
  double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc,CELLDM
  double precision,dimension(3,3), intent(in) :: a
  real(kind(1.d0)), dimension(3) :: xp,xc


  open (25,file='deb_loop.gin_0',form='formatted',status='unknown')
  call lammps_begin ('deb_loop.lmp', nall, a*volc, &
       'Insert_Sias: initial state of the dumbbell jump (deb_loop.gin_0), lengths in Angstrom')


    CELLDM=1.0
  in=1




   write(25,'("1   1   1")')
   write(25,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(25,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(25,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
   write(25,'(i8)') nall

  do i=1,nfix
    ii=nsite_fixed(i)
   do is=1,2
    x = rx(ii)+rsx(is)
    y = ry(ii)+rsy(is)
    z = rz(ii)+rsz(is)

    xp(1:3)=(/x,y,z /)
    xc = MatMul(inv_at(:,:)/volc, xp)*volc
    write(25,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
    call lammps_atom (xp(1), xp(2), xp(3))
   end do
  end do


  do i=1,ntemp
   icard=0
    do is=1,2*nmig
     if (i==nsite_removed(is)) then
      icard=1
     endif
    end do
    do is=1,nfix
     if (i==nsite_fixed(is)) then
      icard=1
     endif
    end do

  if (icard==0) then
    xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
    xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
    write(25,'(a3,3(f18.13),"   1")') '   ' ,xc(1:3)
    call lammps_atom (xp(1), xp(2), xp(3))
   end if

 end do


 icnt=0
! writing sias ....
  do i=1,nmig
   do j=1,3
    xp(1:3)=debatoms(i)%atom(1:3,j)*volc
    xc = MatMul(inv_at(:,:)/volc, xp)*volc
    write(25,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
    call lammps_atom (xp(1), xp(2), xp(3))
   end do
  end do


  call lammps_end ()
  close(25)



  return
  end subroutine deb_write_out_migration

!=======================================================================

   subroutine fin_write_out_migration (nall,ntemp,volc,rx,  &
                  ry,rz,rsx,rsy,rsz,a)
  use math
  use pbc
  use migration_sia
  use the_111_smooth
  use lammps_output

  implicit none
  integer  :: nall,ntemp,i,j,ii,icard,is
  double precision,dimension(nall)     :: rx,ry,rz
  integer  :: in ,icnt
  double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc,CELLDM
  double precision,dimension(3,3), intent(in) :: a
  real(kind(1.d0)), dimension(3) :: xp,xc


  open (26,file='fin_loop.gin_0',form='formatted',status='unknown')
  call lammps_begin ('fin_loop.lmp', nall, a*volc, &
       'Insert_Sias: final state of the dumbbell jump (fin_loop.gin_0), lengths in Angstrom')


    CELLDM=1.0
  in=1



   write(26,'("1   1   1")')
   write(26,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(26,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(26,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
   write(26,'(i8)') nall





  do i=1,nfix
    ii=nsite_fixed(i)
   do is=1,2
    x = rx(ii)+rsx(is)
    y = ry(ii)+rsy(is)
    z = rz(ii)+rsz(is)

    xp(1:3)=(/x,y,z /)
    xc = MatMul(inv_at(:,:)/volc, xp)*volc
    write(26,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
    call lammps_atom (xp(1), xp(2), xp(3))
   end do
  end do

  do i=1,ntemp
   icard=0
    do is=1,2*nmig
     if (i==nsite_removed(is)) then
      icard=1
     endif
    end do
    do is=1,nfix
     if (i==nsite_fixed(is)) then
      icard=1
     endif
    end do

  if (icard==0) then

    xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
    xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
    write(26,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
    call lammps_atom (xp(1), xp(2), xp(3))
   end if
 end do


 icnt=0
! writing sias ....
  do i=1,nmig
   do j=1,3
    xp(1:3)=finatoms(i)%atom(1:3,j)
    xc = MatMul(inv_at(:,:)/volc, xp)
    write(26,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
    call lammps_atom (xp(1), xp(2), xp(3))
   end do
  end do


  call lammps_end ()
  close(26)

  return
  end subroutine fin_write_out_migration

!=======================================================================
