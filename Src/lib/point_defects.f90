!=======================================================================
! point_defects.f90
!
! Point defects: dumbbells, Gao di-interstitials and vacancies. Geometry of
! the inserted atoms and writing of the output files.
!
!   read_sia      vectors of the two dumbbell atoms; template of gao_unit.in
!   write_out     output for dumbbells (types 1, 2, 21, 210, 3, 31)
!   write_out_m   output for Gao di-interstitials (type 4)
!   write_out_vv  output for vacancies (type 5)
!=======================================================================

   subroutine read_sia(tsia,as110,as111,as100,rsx,rsy,rsz,rgx,rgy,rgz)
   use icosaedre, only: rdumb
   use defect_types
   implicit none
   double precision     :: rsx(2),rsy(2),rsz(2),as110,as111,as100
   double precision     :: rgx(6), rgy(6), rgz(6)
   double precision     :: as2, cos01,cos45
   integer :: tsia,i

   cos45 = sqrt(2.d0)/2.d0
   cos01 = sqrt(2.d0)/sqrt(3.d0)

   if (tsia==DEFECT_SIA_110) then
      as2 = as110/2.d0
      rsx(1)= - cos45 * as2
    rsy(1)= - cos45 * as2
    rsz(1)= 0.D0

    rsx(2)= cos45 * as2
    rsy(2)= cos45 * as2
    rsz(2)= 0.d0
   else if ((tsia==DEFECT_SIA_111).or.(tsia==DEFECT_SIA_111_CARE).or. &
            (tsia==DEFECT_SIA_111_MIGRATION).or.(tsia==DEFECT_SIA_111_MIGRATION_CARE)) then
      as2 = as111/2.d0
    rsx(1)= - cos01 * cos45 *  as2
    rsy(1)= - cos01 * cos45 *  as2
    rsz(1)= - cos01 * cos45 *  as2

    rsx(2)= cos01 * cos45 *  as2
    rsy(2)= cos01 * cos45 *  as2
    rsz(2)= cos01 * cos45 *  as2
  else if (tsia==DEFECT_FRANK_LOOP) then

      as2 = as111/2.d0
    rsx(1)= - cos01 * cos45 *  as2 !+ 2.8553/3.0
    rsy(1)= - cos01 * cos45 *  as2 ! -2.8553/6.0
    rsz(1)= - cos01 * cos45 *  as2 ! -2.8553/6.0

    rsx(2)= cos01 * cos45 *  as2 + 2.8553/3.0
    rsy(2)= cos01 * cos45 *  as2  -2.8553/6.0
    rsz(2)= cos01 * cos45 *  as2  -2.8553/6.0

  else if ((tsia==DEFECT_SIA_100).or.(tsia==DEFECT_SIA_100_CARE)) then

      rsx(1)= 0.d0
      rsy(1)= 0.d0
      rsz(1)= -as100/2.d0

      rsx(2)= 0.d0
      rsy(2)= 0.d0
      rsz(2)= as100/2.d0

  else if (tsia==DEFECT_ICO_DUMBBELLS) then
      !rdumb(i,jk,)
      ! i - 1st atom of the dumbell
      ! j - x,y,z of the dumbell
      ! k - the type of the dumbell
      rdumb(1,1,3)= 0.d0
      rdumb(1,2,3)= 0.d0
      rdumb(1,3,3)= -as100/2.d0

      rdumb(2,1,3)= 0.d0
      rdumb(2,2,3)= 0.d0
      rdumb(2,3,3)= as100/2.d0


      rdumb(1,1,2)= 0.d0
      rdumb(1,2,2)= -as100/2.d0
      rdumb(1,3,2)=0.d0

      rdumb(2,1,2)= 0.d0
      rdumb(2,2,2)= as100/2.d0
      rdumb(2,3,2)= 0.d0


      rdumb(1,1,1)= -as100/2.d0
      rdumb(1,2,1)= 0.d0
      rdumb(1,3,1)= 0.d0

      rdumb(2,1,1)= as100/2.d0
      rdumb(2,2,1)= 0.d0
      rdumb(2,3,1)= 0.d0


    else if (tsia==DEFECT_GAO) then

    open (23,file='gao_unit.in',form='formatted' )

    do i=1,6
     read(23,*) rgx(i),rgy(i),rgz(i)
    end do
    close(23)
    else if ((tsia==DEFECT_ICO).or.(tsia==DEFECT_C15).or.(tsia==DEFECT_C15_TWISTED)) then
    write(*,*) 'Bye bye Platon...'
    write(*,*) 'Welcome in the Archimedian solids'


    else
     print*, 'Unknown type of SIA: Check read_sia'
     stop
   end if


   return
   end

!=======================================================================

   subroutine write_out (nall,ntemp,nsia,volc,nsite,rxb,ryb,rzb,rx,  &
                  ry,rz,rsx,rsy,rsz,a,a0)
  use math
  use pbc
  use the_111_smooth
  use lammps_output
  implicit none
  integer  :: nall,ntemp,nsia,n,i,icard,is
  double precision,dimension(nall)     :: rx,ry,rz
  double precision,dimension(ntemp)     :: rxb,ryb,rzb
  integer  :: nsite(nsia),in ,icnt
  double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc,CELLDM
  double precision,dimension(3,3), intent(in) :: a
  real(8) , intent (in) :: a0
  real(8) :: xmin,xmax,ymin,ymax,zmin,zmax,decal
  integer :: inewc
  logical  :: lremove
  real(kind(1.d0)), dimension(3) :: xp,xc
  real(kind(1.d0))               :: xcenter,ycenter,zcenter,adiff, &
                                    xcube,ycube,zcube

  open(11,file='jsia.xyz')
  open(12,file='sia.xyz')
  open(13,file='unit_cell',form='formatted')
  open(15,file='formd_new.str',form='formatted',status='unknown')
  open(16,file='jusia.xyz')
  open(17,file='geom_new.data',form='formatted',status='unknown')
  open(18,file='formd_new_ndm_scl.str',form='formatted',status='unknown')
  open(19,file='usia.xyz')
  open(21,file='grid.gin',form='formatted',status='unknown')
  open(22,file='geom_org.data',form='formatted',status='unknown')
  open(23,file='uvbcsia_sim.xyz')
  open(24,file='uvbcsia_cube.xyz')
  open(26,file='juvsia.xyz')
  call lammps_begin ('formd_new.lmp', nall, a*volc, &
       'Insert_Sias: crystal with defects (formd_new.str), lengths in Angstrom')

    CELLDM=1.0
  in=1



 write(15,'(i8)') nall
   write(15,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(15,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(15,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


   write(17,'(i8)') nall
   write(17,*) 1.0
   write(17,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(17,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(17,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

   write(22,'(i8)') ntemp
   write(22,*) 1.0
   write(22,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(22,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(22,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


   write(18,'("1   1   1")')
   write(18,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(18,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(18,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
   write(18,'(i8)') nall



   write(21,'("1   1   1")')
   write(21,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(21,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(21,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
   write(21,'(i8)') ntemp


10    format(I8,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
     write(13,*) 'number of atoms in the unit cell '
     write(13,*) nall
     write(13,*) 'CELLDM= '
     write(13,*)  CELLDM
   write(13,*) 'position of the atoms in the unit cell'



  write(12,'(i8)') nall
  write(12,'(a)') ' '
  write(11,'(i8)') nall
  write(11,'(a)') ' '
  write(19,'(i8)') 2*nsia
  write(19,'(a)') ' '
  write(16,'(i8)') 2*nsia
  write(16,'(a)') ' '
  write(26,'(i8)') 3*nsia
  write(26,'(a)') ' '
 icnt=0

  xmin=1.d09
  ymin=1.d09
  zmin=1.d09

  xmax=-1.d09
  ymax=-1.d09
  zmax=-1.d09



  do i=1,nsia
       do n=1,2
        icnt =icnt+1
        x = rx(nsite(i))*volc + rsx(n)*volc
        y = ry(nsite(i))*volc + rsy(n)*volc
        z = rz(nsite(i))*volc + rsz(n)*volc
        write(11,'("Fe  ",3f18.13)') x,y,z
        write(12,'("Fe  ",3f18.13)') x,y,z
        write(19,'("Fe  ",3f18.13)') x,y,z
        write(16,'("Fe  ",3f18.13)') x,y,z
        write(26,'("Fe  ",3f18.13)') x,y,z
        write(13,10) icnt,x,y,z,in,in,in
        write(15,'(a3,3(f18.13))') 'Fe ', x,y,z
        call lammps_atom (x, y, z)
        write(17,'(a3,3(f18.13))') '   ', x,y,z
       if (x .ge. xmax) xmax=x
       if (y .ge. ymax) ymax=y
       if (z .ge. zmax) zmax=z
       if (x .le. xmin) xmin=x
       if (y .le. ymin) ymin=y
       if (z .le. zmin) zmin=z
       xp(1:3)=(/  x, y, z/)
       xc = MatMul(inv_at(:,:)/volc, xp)
       write(18,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
       end do
  end do
  !decal=a0/2.d0-0.1d0
  decal=a0+0.1
  xcube=xmax-xmin
  ycube=ymax-ymin
  zcube=zmax-zmin

  xcenter=xmin+xcube/2.d0
  ycenter=ymin+ycube/2.d0
  zcenter=zmin+zcube/2.d0

  adiff=ycube
  if (xcube >= ycube) adiff=xcube
  if (zcube >= adiff)   adiff=zcube


! Write only a slice:
  xmin=xmin-decal
  ymin=ymin-decal
  zmin=zmin-decal
  xmax=xmax+decal
  ymax=ymax+decal
  zmax=zmax+decal


  inewc=0
  do i=1,ntemp
     x=rx(i)*volc
     y=ry(i)*volc
     z=rz(i)*volc
     if (x.gt.xmin) then
      if (y.gt.ymin) then
       if (z.gt.zmin) then
        if (x.lt.xmax) then
         if (y.lt.ymax) then
          if (z.lt.zmax) then
           lremove=.false.
           do is=1,nsia
           if (i==nsite(is))  lremove=.true.
           end do
            if (.not.lremove) inewc=inewc+1
        end if
       end if
      end if
     end if
    end if
   end if
 end do

  write(23,'(i8)') 2*nsia+inewc
  write(23,'(a)') ' '

  do i=1,nsia
       do n=1,2
        icnt =icnt+1
        x = rx(nsite(i))*volc + rsx(n)*volc
        y = ry(nsite(i))*volc + rsy(n)*volc
        z = rz(nsite(i))*volc + rsz(n)*volc
        write(23,'("Fe  ",3f18.13)') x,y,z
       end do
  end do






  do i=1,ntemp
     x=rx(i)*volc
     y=ry(i)*volc
     z=rz(i)*volc
     if (x.gt.xmin) then
      if (y.gt.ymin) then
       if (z.gt.zmin) then
        if (x.lt.xmax) then
         if (y.lt.ymax) then
          if (z.lt.zmax) then
           lremove=.false.
           do is=1,nsia
           if (i==nsite(is))  lremove=.true.
           end do
            if (.not.lremove) write(23,'("Cu  ",3f18.13)') x,y,z
        end if
       end if
      end if
     end if
    end if
   end if
 end do

! Write a cube for Laurent:
  xmin=xcenter-adiff/2.d0-decal
  xmax=xcenter+adiff/2.d0+decal

  ymin=ycenter-adiff/2.d0-decal
  ymax=ycenter+adiff/2.d0+decal

  zmin=zcenter-adiff/2.d0-decal
  zmax=zcenter+adiff/2.d0+decal

  inewc=0
  do i=1,ntemp
     x=rx(i)*volc
     y=ry(i)*volc
     z=rz(i)*volc
     if (x.gt.xmin) then
      if (y.gt.ymin) then
       if (z.gt.zmin) then
        if (x.lt.xmax) then
         if (y.lt.ymax) then
          if (z.lt.zmax) then
           lremove=.false.
           do is=1,nsia
           if (i==nsite(is))  lremove=.true.
           end do
            if (.not.lremove) inewc=inewc+1
        end if
       end if
      end if
     end if
    end if
   end if
 end do

  write(24,'(i8)') 2*nsia+inewc
  write(24,'(a)') ' '

  do i=1,nsia
       do n=1,2
        icnt =icnt+1
        x = rx(nsite(i))*volc + rsx(n)*volc
        y = ry(nsite(i))*volc + rsy(n)*volc
        z = rz(nsite(i))*volc + rsz(n)*volc
        write(24,'("Fe  ",3f18.13)') x,y,z
       end do
  end do






  do i=1,ntemp
     x=rx(i)*volc
     y=ry(i)*volc
     z=rz(i)*volc
     if (x.gt.xmin) then
      if (y.gt.ymin) then
       if (z.gt.zmin) then
        if (x.lt.xmax) then
         if (y.lt.ymax) then
          if (z.lt.zmax) then
           lremove=.false.
           do is=1,nsia
           if (i==nsite(is))  lremove=.true.
           end do
            if (.not.lremove) write(24,'("Cu  ",3f18.13)') x,y,z
        end if
       end if
      end if
     end if
    end if
   end if
 end do






  do i=1,ntemp
   xp(1:3)=(/ rxb(i),ryb(i),rzb(i) /)*volc
   xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
   write(21,'(a3,3(f18.13),"   1")') '   ' ,xc(1:3)
   write(22,'(a3,3(f18.13),"   1")')  '   ',rxb(i)*volc,ryb(i)*volc,rzb(i)*volc
   xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
   xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))

   icard=0
   do is=1,nsia
     if (i==nsite(is)) then
      icard=1
     endif
   end do











  if (icard==0) then
   icnt=icnt+1
   write(11,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(12,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(13,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
   write(15,'(a3,3(f18.13))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
   call lammps_atom (rx(i)*volc, ry(i)*volc, rz(i)*volc)
   write(17,'(a3,3(f18.13))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(18,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
  end if
  if (icard==1)  write(26,'(a3,3(f18.13),3i6)') 'He ', &
                 rx(i)*volc,ry(i)*volc,rz(i)*volc, &
                 nint(rx(i)/(a0/2.d0)), &
                 nint(ry(i)/(a0/2.d0)), &
                 nint(rz(i)/(a0/2.d0))
 end do


  open(14,file='unit_vect',form='formatted')
  write (14,*) ' number of unit vectors 0 1 2 ou 3'
  write(14,*) 3
  write(14,*) 'unit vectors:'
  write(14,'(3f18.13)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
  write(14,'(3f18.13)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
  write(14,'(3f18.13)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

  write(19,*) ' '
  write(19,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
  write(19,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
  write(19,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

  write(12,*) ' '
  write(12,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
  write(12,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
  write(12,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc




  call lammps_end ()

  close(11)
  close(12)
  close(13)
  close(14)
  close(16)
  close(19)
  close(23)
  close(24)




  return
  end

!=======================================================================

   subroutine write_out_m (nall,ntemp,nsia,volc,a0,nsite_m,rx,  &
                   ry,rz,rsx,rsy,rsz,a)
   use math
  use pbc
   use lammps_output
   integer  :: nall,ntemp,nsia,n,i,icard,in,icnt,is,ii
   double precision     :: rx(nall),ry(nall),rz(nall),a0,a02,CELLDM
   integer  :: nsite_m(4,nsia)
   double precision     :: rsx(6),rsy(6),rsz(6),x,y,z, volc
   double precision,dimension(3,3) :: a
   real(kind(1.d0)),dimension(3)  :: xp,xc

   open(11,file='jsia.xyz')
   open(12,file='sia.xyz')
   open(13,file='unit_cell',form='formatted')
   open(15,file='formd_new.str',form='formatted',status='unknown')
   open(16,file='jusia.xyz')
   open(17,file='geom_new.data',form='formatted',status='unknown')
   open(18,file='formd_new_ndm_scl.str',form='formatted',status='unknown')
   open(19,file='usia.xyz')
   open(21,file='formd_new_ndm_org.str',form='formatted',status='unknown')
   open(22,file='geom_org.data',form='formatted',status='unknown')
   call lammps_begin ('formd_new.lmp', nall, a*volc, &
        'Insert_Sias: crystal with defects (formd_new.str), lengths in Angstrom')

     CELLDM=1.0
   in=1

   a02=a0/2.d0


  write(15,'(i8)') nall
    write(15,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(15,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(15,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


  write(17,'(i8)') nall
  write(17,*) 1.0
    write(17,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(17,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(17,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

  write(22,'(i8)') ntemp
  write(22,*) 1.0
    write(22,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(22,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(22,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


    write(18,'("1   1   1")')
  write(18,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(18,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(18,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
  write(18,'(i8)') nall

    write(21,'("1   1   1")')
  write(21,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(21,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(21,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
  write(21,'(i8)') ntemp


10       format(I8,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
      write(13,*) 'number of atoms in the unit cell '
      write(13,*) nall
      write(13,*) 'CELLDM= '
      write(13,*)  CELLDM
    write(13,*) 'position of the atoms in the unit cell'



   write(12,'(i8)') nall
   write(12,'(a)') ' '
   write(11,'(i8)') nall
   write(11,'(a)') ' '
   write(19,'(i8)') 6*nsia
   write(19,'(a)') ' '
   write(16,'(i8)') 6*nsia
   write(16,'(a)') ' '
   icnt=0
   do i=1,nsia
    do n=1,6
   icnt =icnt+1
     x = rx(nsite_m(1,i))*volc + rsx(n)*volc*a0
     y = ry(nsite_m(1,i))*volc + rsy(n)*volc*a0
     z = rz(nsite_m(1,i))*volc + rsz(n)*volc*a0
   write(11,'("Fe  ",3f18.13)') x,y,z
   write(12,'("Fe  ",3f18.13)') x,y,z
   write(19,'("Fe  ",3f18.13)') x,y,z
   write(16,'("Fe  ",3f18.13)') x,y,z
   write(13,10) icnt,x,y,z,in,in,in
   write(15,'(a3,3(f18.13))') 'Fe ', x,y,z
   call lammps_atom (x, y, z)
   write(17,'(a3,3(f18.13))') '   ', x,y,z
    xp(1:3)=(/  x, y, z/)
    xc(1:3) = MatMul(inv_at(1:3,1:3)/volc, xp(1:3))
    write(18,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
  end do
   end do


   do i=1,ntemp
    xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
    xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
    write(21,'(a3,3(f18.13),"   1")') '   ' ,xc(1:3)

    write(22,'(a3,3(f18.13),"   1")')  '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
    icard=0
    do is=1,nsia
     do ii=1,4
      if (i==nsite_m(ii,is)) then
       icard=1
      endif
     end do
  end do
   if (icard==0) then
   icnt=icnt+1
   write(11,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(12,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(13,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
   write(15,'(a3,3(f18.13))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
   call lammps_atom (rx(i)*volc, ry(i)*volc, rz(i)*volc)
   write(17,'(a3,3(f18.13))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(18,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
   end if
   end do


   open(14,file='unit_vect',form='formatted')
   write (14,*) ' number of unit vectors 0 1 2 ou 3'
   write(14,*) 3
   write(14,*) 'unit vectors:'
   write(14,'(3f18.13)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(14,'(3f18.13)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(14,'(3f18.13)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

   write(19,*) ' '
   write(19,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(19,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(19,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

   write(12,*) ' '
   write(12,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(12,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(12,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc




   call lammps_end ()

   close(12)
   close(11)
   close(13)
   close(14)
   close(16)
   close(19)



   return
   end

!=======================================================================

   subroutine write_out_vv (nall,ntemp,nsia,nsize,volc,nsite,rx,  &
                   ry,rz,a)
   use math
   use pbc
   use lammps_output
   integer  :: nall,ntemp,nsia,i,icard,nsize,in,icnt,is
   double precision     :: rx(nsize),ry(nsize),rz(nsize),CELLDM
   integer  :: nsite(nsia)
   double precision     :: x,y,z, volc
   double precision,dimension(:,:) :: a(3,3)
   real(kind(1.d0)),dimension(3) :: xp,xc

   open(11,file='jsia.xyz')
   open(12,file='sia.xyz')
   open(19,file='usia.xyz')
   open(16,file='jusia.xyz')
   open(16,file='jusia.xyz')
   open(13,file='unit_cell',form='formatted')
   open (17,file='geom_new.data',form='formatted',status='unknown')
   open (22,file='geom_org.data',form='formatted',status='unknown')
   open (18,file='formd_new_ndm_scl.str',form='formatted',status='unknown')
   open (21,file='grid.gin',form='formatted',status='unknown')

     CELLDM=1.0
   in=1


  open (15,file='formd_new.str',form='formatted',status='unknown')
  call lammps_begin ('formd_new.lmp', nall, a*volc, &
       'Insert_Sias: crystal with defects (formd_new.str), lengths in Angstrom')

  write(15,'(i8)') nall
    write(15,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(15,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(15,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


  write(17,'(i8)') nall
  write(17,*) 1.0
    write(17,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(17,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(17,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc



    write(18,'("1   1   1")')
  write(18,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(18,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(18,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
  write(18,'(i8)') nall

  write(22,'(i8)') ntemp
  write(22,*) 1.0
    write(22,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(22,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(22,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

    write(21,'("1   1   1")')
  write(21,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
    write(21,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
    write(21,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
  write(21,'(i8)') ntemp


10       format(I8,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
      write(13,*) 'number of atoms in the unit cell '
      write(13,*) nall
      write(13,*) 'CELLDM= '
      write(13,*)  CELLDM
    write(13,*) 'position of the atoms in the unit cell'



   write(11,'(i8)') ntemp
   write(12,'(i8)') ntemp
   write(11,'(a)') ' '
   write(12,'(a)') ' '
   write(19,'(i8)') nsia
   write(19,'(a)') ' '
   write(16,'(i8)') nsia
   write(16,'(a)') ' '
   icnt=0
   do i=1,nsia
     x = rx(nsite(i))*volc
     y = ry(nsite(i))*volc
     z = rz(nsite(i))*volc
    write(11,'("Fe  ",3f18.13)') x,y,z
    write(12,'("Fe  ",3f18.13)') x,y,z
    write(19,'("Fe  ",3f18.13)') x,y,z
    write(16,'("Fe  ",3f18.13)') x,y,z
   end do


   do i=1,ntemp
    xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
    xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
    write(21,'(a3,3(f18.13),"   1")') '   ' ,xc(1:3)
   write(22,'(a3,3(f18.13),"   1")')  '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
   icard=0
    do is=1,nsia
     if (i==nsite(is)) then
     icard=1
     endif
  end do
   if (icard==0) then
   icnt=icnt+1
   write(11,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(12,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(13,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
   write(15,'(a3,3(f18.13))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
   call lammps_atom (rx(i)*volc, ry(i)*volc, rz(i)*volc)
   write(17,'(a3,3(f18.13))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
   write(18,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
   end if
   end do

   write(*,*) 'Final countdown',  icnt,ntemp,nall
   open(14,file='unit_vect',form='formatted')
   write (14,*) ' number of unit vectors 0 1 2 ou 3'
   write(14,*) 3
   write(14,*) 'unit vectors:'
   write(14,'(3f18.13)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(14,'(3f18.13)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(14,'(3f18.13)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

   write(19,*) ' '
   write(19,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(19,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(19,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

   write(12,*) ' '
   write(12,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(12,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(12,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc




   call lammps_end ()

   close(12)
   close(11)
   close(13)
   close(14)
   close(16)
   close(19)



   return
   end

!=======================================================================
