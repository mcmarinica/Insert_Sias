!#####################################################################
   
   subroutine read_ntemp (ntemp)
    integer :: ntemp
     
   open (10,file='formd.str',form='formatted',status='unknown')
     read (10,'(i8)') ntemp
    
    return
    end
!#####################################################################
  
   subroutine read_nsia (nsia,tsia)
   use migration_sia
   integer :: nsia,tsia
    
   if ((tsia==22).or.(tsia==221)) then
    open (11,file='no00m.inp',form='formatted',status='unknown')
    read (11,*) nsia,nmig,nfix
    if ((nfix+nmig)/=nsia) then
      write(*,*) "nfix + nmig should be equal to nsia, nsia,nfix,nmig.....:",nsia,nfix,nmig
      write(*,*) "stop in read_nsia"
      stop
    end if  
   
   else 
  
    open (11,file='no00.inp',form='formatted',status='unknown')
    read (11,*) nsia
  
   end if 
  
   return
   end
           
!#####################################################################
   subroutine read_input(ntemp,nsize,a,rx,ry,rz)      
   use math
   use pbc
   implicit none
   integer              :: ntemp,nsize,i
   double precision     :: rx(nsize),ry(nsize),rz(nsize)
   double precision,dimension(:,:) :: a(3,3)
   character*3 :: sym
   

    read(10,'(3f18.13)')a(1,1), a(1,2), a(1,3)
    read(10,'(3f18.13)')a(2,1), a(2,2), a(2,3)
    read(10,'(3f18.13)')a(3,1), a(3,2), a(3,3)

    at(:,1)=a(1,:)
    at(:,2)=a(2,:)
    at(:,3)=a(3,:)
    ! Inverse matrix
    CALL MatInv(at,inv_at)
  
  do i=1,ntemp
        read(10,'(a3,3(f18.13))') sym, rx(i),ry(i),rz(i)
    end do             
 
        
   
   return
   end 
!#####################################################################       
   subroutine read_config_sia (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
   integer  :: nsia, n1MAX,n2MAX,n3MAX,i  
   integer  :: ns1(nsia),ns2(nsia),ns3(nsia)
        
      
      do i=1,nsia
       read (11,*) ns1(i),ns2(i),ns3(i)
       !write(*,*) 'ns11',ns1(i),ns2(i),ns3(i)
        if ( (ns1(i) > n1MAX+1) .or.(ns2(i) > n2MAX+1) .or. &
          (ns3(i) > n3MAX+1) ) then
        print*, 'Big problem in read_config_sia'
        print*, 'The cell is not big enough. Increase the size!'
        print*, '<read_config_sia>', n1MAX+1, n2MAX+1, n3MAX+1
        stop
        end if          
      end do 
   return
   end 

!#####################################################################       
   subroutine read_config_sia_migration (nsia, ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
   use migration_sia 
   integer  :: nsia, n1MAX,n2MAX,n3MAX,i  
   integer  :: ns1(nsia),ns2(nsia),ns3(nsia)
   integer  :: i1,i2     
     i1=0
     i2=0

      do i=1,nsia
       read (11,*) ns1(i),ns2(i),ns3(i),type_migration(i)
       !write(*,*) 'ns11',ns1(i),ns2(i),ns3(i)
        if ( (ns1(i) > n1MAX+1) .or.(ns2(i) > n2MAX+1) .or. &
          (ns3(i) > n3MAX+1) ) then
        print*, 'Big problem in read_config_sia'
        print*, 'The cell is not big enough. Increase the size!'
        print*, '<read_config_sia_migration>'
        stop
        end if
        if (.not.((type_migration(i)==0).or.(type_migration(i)==1))) then
          write(*,*) "The 4th colomn of no00m.inp should be 0 or 1"
          write(*,*) "stop in read_config_sia_migration"
          stop
        end if

        if (type_migration(i)==1) then
         if ( (ns1(i) > n1MAX-1) .or.(ns2(i) > n2MAX-1) .or. &
           (ns3(i) > n3MAX-1) ) then
           print*, 'Big problem in read_config_sia_migration'
           print*, 'The cell is not big enough. Increase the size!'
           print*, '<read_config_sia_migration>'
           stop
         end if
         if ( (ns1(i)==0) .or.(ns2(i)==0) .or. &
           (ns3(i)==0 )) then
           print*, 'Big problem in read_config_sia_migrration'
           print*, 'You cannot have 0 in the coordinates for migrating interstitials'
           print*, '<read_config_sia_migration>'
           stop
         end if
        end if   
       
                  
       if (type_migration(i)==1) i1=i1+1
       if (type_migration(i)==0) i2=i2+1
      end do 
            
     if (i1/=nmig) then
      write(*,*) "The number of migrating atoms diff from nmig ...i1 and nmig ...:", i1,nmig
      stop
     end if

     if (i2/=nfix) then
      write(*,*) "The number of fixed atoms diff from nfix ...i2 and nfix ...:", i2,nfix
      stop
     end if



   return
   end 
           
          
!##################################################################### 

      
   subroutine analyze_config_sia (tsia,nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
   
   implicit none 
   integer  :: nsia, n1MAX,n2MAX,n3MAX,i,j  
   integer  :: ns1(nsia),ns2(nsia),ns3(nsia)
   integer  :: ii,neigh(5),tsia
   real(8)  :: bcc(5),dij,delta=0.1d0 

  if ((tsia>3).and.(tsia/=21)) then
   write(*,*) 'WARNING ... the analyze_config SIA .... Not yet implemented for tsia  ',tsia
   !stop
  end if  

  if (tsia==2) then
 ! in a0 units the bcc NN distances
  bcc(1)=dsqrt(3.d0)/2.d0
  bcc(2)=1.d0
  bcc(3)=dsqrt(2.d0)
  bcc(4)=dsqrt(3.d0)
  bcc(5)=2.d0
 end if 

if (tsia==3) then
  bcc(1)=dsqrt(3.d0)/2.d0
  bcc(2)=1.d0
  bcc(3)=dsqrt(2.d0)
  bcc(4)=dsqrt(11.d0)/2.d0
  bcc(5)=2.d0
end if 

 ! in internal units of a0/2 this distances gives
  bcc(:)=2.d0*bcc(:)
       
 write(567,'(i8)')  nsia
 do i=1,nsia
    neigh(:) = 0
  do j=1,nsia
  ! if (tsia==3) then
  !  dij=dsqrt(dble(                      &
  !            (ns1(i)-ns1(j))**2 +  &
  !            (ns2(i)-ns2(j))**2    &
  !            ))
  !
  !else 
    dij=dsqrt(dble(                      &
              (ns1(i)-ns1(j))**2 +  &
              (ns2(i)-ns2(j))**2 +  &
              (ns3(i)-ns3(j))**2    &
              ))
   !end if
    do ii=1,5
     if ((dij<=(bcc(ii)+delta)).and.(dij>=(bcc(ii)-delta))) then
        neigh(ii)=neigh(ii)+1
     end if
    end do
   end do !j-nsia
   write(567,'(i8,5i6)') i,neigh(1:5) 
 end do  !i-nsia
            
   return
   end 
           
!#####################################################################      
   subroutine write_n123(ntemp,nsize,a0,n1,n2,n3,rx,ry,rz,  &
                         n1MAX,n2MAX,n3MAX)      
   use math
   use the_111_smooth, ONLY : pbc_111
   use the_100_smooth, ONLY : pbc_100
   implicit none

   integer :: nsize, i,n1MAX,n2MAX,n3MAX,ntemp
   integer  :: n1(nsize),n2(nsize),n3(nsize)
   real(8)     :: rx(nsize),ry(nsize),rz(nsize)
   integer :: itemp,ngcd,n1pbc,n2pbc,n3pbc
   
   double precision     :: a0,a02
   
   a02= a0/2.d0
   
     do i=1,ntemp
          n1(i)=nint(rx(i)/a02)
          n2(i)=nint(ry(i)/a02)
          n3(i)=nint(rz(i)/a02)
     end do
    n1MAX = maxval(n1)
    n2MAX=  maxval(n2) 
    n3MAX=  maxval(n3) 

    n1pbc=n1MAX+1
    n2pbc=n2MAX+1
    n3pbc=n3MAX+1
    
    print*, 'n1MAX,n2MAX,n3MAX',n1MAX,n2MAX,n3MAX

    call gcd (n1pbc,n2pbc,itemp)
    call gcd (itemp,n3pbc,ngcd)
    pbc_111=n1pbc*(n2pbc/ngcd)*(n3pbc/ngcd)
    pbc_100= n3pbc
     write(*,*) 'We will have periodicity in <111> direction in', pbc_111
     write(*,*) 'We will have periodicity in <100> direction in', pbc_100
  
   return
   end 
!#####################################################################



   subroutine read_sia(tsia,as110,as111,as100,rsx,rsy,rsz,rgx,rgy,rgz)      
   implicit none
   double precision     :: rsx(2),rsy(2),rsz(2),as110,as111,as100
   double precision     :: rgx(6), rgy(6), rgz(6)
   double precision     :: as2, cos01,cos45
   integer :: tsia,i
   
   cos45 = sqrt(2.d0)/2.d0
   cos01 = sqrt(2.d0)/sqrt(3.d0)
   
   if (tsia==1) then 
      as2 = as110/2.d0
      rsx(1)= - cos45 * as2
    rsy(1)= - cos45 * as2
    rsz(1)= 0.D0
    
    rsx(2)= cos45 * as2
    rsy(2)= cos45 * as2
    rsz(2)= 0.d0
   else if ((tsia==2).or.(tsia==22).or.(tsia==21).or.(tsia==221)) then
      as2 = as111/2.d0
    rsx(1)= - cos01 * cos45 *  as2
    rsy(1)= - cos01 * cos45 *  as2
    rsz(1)= - cos01 * cos45 *  as2
    
    rsx(2)= cos01 * cos45 *  as2
    rsy(2)= cos01 * cos45 *  as2
    rsz(2)= cos01 * cos45 *  as2


    
   else if ((tsia==3).or.(tsia==31)) then
      rsx(1)= 0.d0
    rsy(1)= 0.d0
    rsz(1)= -as100/2.d0
    
    rsx(2)= 0.d0
    rsy(2)= 0.d0
    rsz(2)= as100/2.d0
    else if (tsia==4) then
    
    open (23,file='gao_unit.in',form='formatted' )
     
    do i=1,6 
     read(23,*) rgx(i),rgy(i),rgz(i)
    end do
    close(23)
    else if ((tsia==7).or.(tsia==6).or.(tsia==67)) then
    write(*,*) 'Bye bye Platon...'
    write(*,*) 'Welcome in the Archimedian solids'
   

    else
     print*, 'Unknown type of SIA: Check read_sia' 
     stop
   end if    
        
    
   return
   end 
!#####################################################################
   subroutine write_out (nall,ntemp,nsia,volc,nsite,rxb,ryb,rzb,rx,  &
                  ry,rz,rsx,rsy,rsz,a,a0)
  use math
  use pbc 
  use the_111_smooth
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
  write(*,*) xmin, xmax
  write(*,*) xcenter,ycenter,zcenter
  write(*,*) adiff 


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
!#####################################################################
       
!#####################################################################
               
   subroutine write_out_m (nall,ntemp,nsia,volc,a0,nsite_m,rx,  &
                   ry,rz,rsx,rsy,rsz,a)
   use math
  use pbc
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

  


   close(12)
   close(11)
   close(13)
   close(14)
   close(16)
   close(19)
   


   return
   end
!#####################################################################
           
   subroutine write_out_vv (nall,ntemp,nsia,nsize,volc,nsite,rx,  &
                   ry,rz,a)
   use math
   use pbc
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

  


   close(12)
   close(11)
   close(13)
   close(14)
   close(16)
   close(19)
   


   return
   end
!#####################################################################



!#####################################################################
   subroutine write_out_migration (nall,ntemp,volc,rx,  &
                  ry,rz,rsx,rsy,rsz,a)
  use math
  use pbc
  use migration_sia 
  use the_111_smooth

  implicit none
  integer  :: nall,ntemp,i,j,ii,icard,is
  double precision,dimension(nall)     :: rx,ry,rz
  integer  :: in ,icnt      
  double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc,CELLDM       
  double precision,dimension(3,3), intent(in) :: a
  real(kind(1.d0)), dimension(3) :: xp,xc

  
  open (25,file='deb_loop.gin_0',form='formatted',status='unknown')
  open (26,file='fin_loop.gin_0',form='formatted',status='unknown')


    CELLDM=1.0
  in=1
  



   write(25,'("1   1   1")') 
   write(25,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(25,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(25,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
   write(25,'(i8)') nall

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
    write(25,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
    write(26,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
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
    write(26,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
   end if

 end do


 icnt=0
! writing sias ....
  do i=1,nmig
   do j=1,3
    xp(1:3)=debatoms(i)%atom(1:3,j)*volc
    xc = MatMul(inv_at(:,:)/volc, xp)*volc
    write(25,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
   end do
   do j=1,3
    xp(1:3)=finatoms(i)%atom(1:3,j)
    xc = MatMul(inv_at(:,:)/volc, xp)
    write(26,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
   end do
  end do
   
  
  close(25)
  close(26)
  


  return
  end
!#####################################################################
   
!#####################################################################
   subroutine deb_write_out_migration (nall,ntemp,volc,rx,  &
                  ry,rz,rsx,rsy,rsz,a)
  use math
  use pbc
  use migration_sia 
  use the_111_smooth

  implicit none
  integer  :: nall,ntemp,i,j,ii,icard,is
  double precision,dimension(nall)     :: rx,ry,rz
  integer  :: in ,icnt      
  double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc,CELLDM       
  double precision,dimension(3,3), intent(in) :: a
  real(kind(1.d0)), dimension(3) :: xp,xc

  
  open (25,file='deb_loop.gin_0',form='formatted',status='unknown')


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
    if (i==7) write(*,*) "i",  xc(1:3)
   end if

 end do


 icnt=0
! writing sias ....
  do i=1,nmig
   do j=1,3
    xp(1:3)=debatoms(i)%atom(1:3,j)*volc
    xc = MatMul(inv_at(:,:)/volc, xp)*volc
    write(25,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
   end do
  end do
   
  
  close(25)
  


  return
  end subroutine deb_write_out_migration 
!#####################################################################
!#####################################################################
   subroutine fin_write_out_migration (nall,ntemp,volc,rx,  &
                  ry,rz,rsx,rsy,rsz,a)
  use math
  use pbc
  use migration_sia 
  use the_111_smooth

  implicit none
  integer  :: nall,ntemp,i,j,ii,icard,is
  double precision,dimension(nall)     :: rx,ry,rz
  integer  :: in ,icnt      
  double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc,CELLDM       
  double precision,dimension(3,3), intent(in) :: a
  real(kind(1.d0)), dimension(3) :: xp,xc

  
  open (26,file='fin_loop.gin_0',form='formatted',status='unknown')


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
    if (i==7) write(*,*) "f",  xc(1:3)
   end if
 end do


 icnt=0
! writing sias ....
  do i=1,nmig
   do j=1,3
    xp(1:3)=finatoms(i)%atom(1:3,j)
    xc = MatMul(inv_at(:,:)/volc, xp)
    write(26,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
   end do
  end do
   
  
  close(26)

  return
  end subroutine fin_write_out_migration
!#####################################################################

