!#####################################################################
               
       
       subroutine read_config_sia_octa (nsia,type_octa,  &
                         ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
       
       integer  :: nsia, n1MAX,n2MAX,n3MAX,i  
       integer  :: ns1(nsia),ns2(nsia),ns3(nsia),type_octa(nsia)
            
          
          do i=1,nsia
           read (11,*) ns1(i),ns2(i),ns3(i)
           !write(*,*) 'ns11',ns1(i),ns2(i),ns3(i)
            if ( (ns1(i) > n1MAX) .or.(ns2(i) > n2MAX) .or. &
              (ns3(i) > n3MAX) ) then
            print*, 'Big problem in read_config_sia'
            print*, 'The cell is not big enough. Increase the size!'
              write(*,*) '<read_config_sia_octa>'
            stop
            end if          
          end do 
          do i=1,nsia
           read(11,*) type_octa(i)
           if ( (type_octa(i).lt.1).or.((type_octa(i).gt.3)) ) then
           write(*,*) 'There is no implementation for this octa type', &
             type_octa(i)
            write(*,*) '<read_config_sia_octa>'
            stop
           end if
          end do
                
       return
       end 
!#####################################################################

      subroutine  lattice_compute_different_removeatoms (nsia, & 
          nlac_remove, nsite_octa_vac_remove,indx,list_no_order,ndatoms)
      
      
      implicit none
        !intent IN/OUT
      integer, intent(in) :: nsia,nlac_remove
        integer, intent(out):: ndatoms
        integer, intent(out) :: indx(nlac_remove)
      real(8), intent(out)  :: list_no_order(nlac_remove)
      !local variables
        integer :: ii,i_double,icnt,i
      integer :: nsite_octa_vac_remove(10*nsia)
      
      if (nlac_remove.gt.10*nsia) then
       write(*,*) 'nlac_remove is greater than 10*nsia'
       write(*,*) 'Problems!!'
       stop
      end if
      
        icnt=0
        do ii=1,nlac_remove
          list_no_order(ii)=nsite_octa_vac_remove(ii)
        end do
      
       call indexx (nlac_remove,list_no_order,indx) 
       

       i_double=0
       !debug write(*,*) 'list', list_no_order(indx(1)),i_double
       do    i=1,nlac_remove-1
         if (dabs(list_no_order(indx(i+1))-list_no_order(indx(i))).lt.0.01d0) then 
        i_double = i_double+1 
       end if
       !debug write(*,*) 'list', list_no_order(indx(i+1)),i_double
       end do
       
        do ii=1,nlac_remove
         if(dabs(list_no_order(ii)-nsite_octa_vac_remove(ii)).gt.1d-5)& 
             then
           write(*,*) '!!!! in list_no_order and nsite_octa_vac_remove '
            end if
        end do

       
         write(*,*) 'The number of vacancies which are superposed....... ........: ',nlac_remove        
         write(*,*) 'The number of DIFFERENT VACANCIES left in the input ........: ',nlac_remove-i_double
         ndatoms=nlac_remove-i_double
       
       return
       end  
         
!#####################################################################
               
      subroutine  lattice_compute_atoms_to_remove  (atoms_to_remove,&
                    atoms_to_remove_type,ndatoms,nlac_remove,         &
                     list_no_order,nsite_octa_vac_remove_type,indx)

         integer, intent(in)  :: ndatoms,nlac_remove
       real(8), intent(in)  :: list_no_order(nlac_remove)
       integer, intent(in)  :: indx(nlac_remove),nsite_octa_vac_remove_type(nlac_remove)
       integer, intent(out)  :: atoms_to_remove_type(ndatoms)
       real(8), intent(out)   :: atoms_to_remove(ndatoms)
         integer :: i_double,icnt,i 
               
       
       atoms_to_remove(1)=list_no_order(indx(1))
       atoms_to_remove_type(1)=nsite_octa_vac_remove_type(indx(1))
       i_double=0
       icnt=1
       do    i=1,nlac_remove-1
         if (dabs(list_no_order(indx(i+1))-list_no_order(indx(i))).lt.0.01d0) then 
        i_double = i_double+1 
       else 
        icnt=icnt+1
        atoms_to_remove(icnt)=list_no_order(indx(i+1))
        atoms_to_remove_type(icnt)=nsite_octa_vac_remove_type(indx(i+1))
        !debug write(*,*) 'toremove', atoms_to_remove(icnt)
       end if 
       end do
       !debug write(*,*) 'list', icnt,i_double,10*nsia-i_double
       
       
       return
       end
!#####################################################################

      subroutine  offlattice_get_coordinates_and_order(a0,rxgaos,rygaos,rzgaos,&
                                                       rx,    ry,    rz,    &
                     offlattice_no_order,indx,nsite,type_octa,nsia,nsize,gndatoms,natoms_extra)

       use octahedre
       implicit none

       integer, intent(in) :: nsia,nsize
       integer, intent(in) :: nsite(nsia),type_octa(nsia)
       real(8), intent(in)   :: a0
       real(8), dimension(nsize)  :: rx,ry,rz
       integer,  intent(out) :: indx (12*nsia),gndatoms,natoms_extra
       
       real(8), dimension(12*nsia),  intent(out) :: offlattice_no_order,rxgaos,rygaos,rzgaos
       integer :: i,it_octo,i_double,icnt,ii,jj
       real(8)   :: a02,alpha,origin(3)
               
             
       a02=a0/2.d0
       alpha=0.5d0
       icnt=0
      do i=1,nsia
       it_octo=type_octa(i)
         origin(1:3)=(/ rx(nsite(i)),  ry(nsite(i)), rz(nsite(i)) /)   
           if ( (it_octo.eq.1).or.(it_octo.eq.2) )then
          do ii=1,4  !loop over the 4 corners of the tetrahedron
           do jj=1,3 !loop over the three atoms located in the corner of the tetrahedron
               icnt=icnt+1
             rxgaos(icnt)=a02*(alpha*octo(it_octo)%gaos(1,jj,ii)+octo(it_octo)%atom(1,ii))+ origin(1)
             rygaos(icnt)=a02*(alpha*octo(it_octo)%gaos(2,jj,ii)+octo(it_octo)%atom(2,ii))+ origin(2)
             rzgaos(icnt)=a02*(alpha*octo(it_octo)%gaos(3,jj,ii)+octo(it_octo)%atom(3,ii))+ origin(3)
             offlattice_no_order(icnt)=rxgaos(icnt)**3+10*rygaos(icnt)**3+100*rzgaos(icnt)**3
           end do
          end do
          else if (it_octo.eq.4) then
          
          do ii=1,1  !loop over the 4 corners of the tetrahedron
           do jj=1,3 !loop over the three atoms located in the corner of the tetrahedron
               icnt=icnt+1
             rxgaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(1,jj,ii)+octo(it_octo-2)%atom(1,ii))+ origin(1)
             rygaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(2,jj,ii)+octo(it_octo-2)%atom(2,ii))+ origin(2)
             rzgaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(3,jj,ii)+octo(it_octo-2)%atom(3,ii))+ origin(3)
             offlattice_no_order(icnt)=rxgaos(icnt)**3+10*rygaos(icnt)**3+100*rzgaos(icnt)**3
           end do
          end do
          
          do ii=3,3  !loop over the 4 corners of the tetrahedron
           do jj=2,2 !loop over the three atoms located in the corner of the tetrahedron
               icnt=icnt+1
             rxgaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(1,jj,ii)+octo(it_octo-2)%atom(1,ii))+ origin(1)
             rygaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(2,jj,ii)+octo(it_octo-2)%atom(2,ii))+ origin(2)
             rzgaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(3,jj,ii)+octo(it_octo-2)%atom(3,ii))+ origin(3)
             offlattice_no_order(icnt)=rxgaos(icnt)**3+10*rygaos(icnt)**3+100*rzgaos(icnt)**3
           end do
          end do

          end if    
      end do 
         
       natoms_extra=icnt
!old       call indexx (12*nsia, offlattice_no_order,indx) 
       call indexx (icnt, offlattice_no_order,indx) 

       i_double=0
       !debug write(*,'("listgaos",f14.3,i4,3f12.3)') &
       !debug offlattice_no_order(indx(1)),i_double,rxgaos(1),rygaos(1),rzgaos(1)
       do    i=1,icnt-1
         !old was 1.1
       if (dabs(offlattice_no_order(indx(i+1))-offlattice_no_order(indx(i))).lt.1.1d0) then 
       i_double = i_double+1 
       end if
       !debug write(*,'("listgaos",f14.3,i4,3f12.3)') &
       !debug offlattice_no_order(indx(i+1)),i_double,&
       !debug rxgaos(indx(i+1)),rygaos(indx(i+1)),rzgaos(indx(i+1))
       end do
       
       
         write(*,*) 'The number of DIFFERENT GAOS atoms which must be removed ....: ',i_double        
         write(*,*) 'The number of DIFFERENT GAOS atoms left in input ............: ',icnt-i_double
         gndatoms=icnt-i_double

       return
       end

       subroutine offlattice_compute_atoms_to_remove(gndatoms,rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos,&
                             offlattice_no_order,indx,nsia,natoms_extra)

       implicit none

       integer, intent(in) :: nsia
       integer,  intent(in) :: indx (12*nsia),gndatoms,natoms_extra
       
       real(8), dimension(12*nsia),    intent(in) :: offlattice_no_order,temprxgaos,temprygaos,temprzgaos
       real(8), dimension(gndatoms),  intent(out) :: rxgaos,rygaos,rzgaos
       integer :: i_double,i,icnt

       i_double=0
       icnt=1
       
        rxgaos(1)=temprxgaos(indx(1))
        rygaos(1)=temprygaos(indx(1))
        rzgaos(1)=temprzgaos(indx(1))
       do    i=1,natoms_extra-1
         if (dabs(offlattice_no_order(indx(i+1))-offlattice_no_order(indx(i))).lt.0.1d0) then 
       i_double = i_double+1 
       else
        icnt=icnt+1
        rxgaos(icnt)=temprxgaos(indx(i+1))
        rygaos(icnt)=temprygaos(indx(i+1))
        rzgaos(icnt)=temprzgaos(indx(i+1))

       end if
       end do
       
         write(*,'("II The number of DIFFERENT GAOS atoms which &
           must be removed .: ",i6)') ,i_double        
       return
       end
       
!#####################################################################
               



       subroutine write_coordinates (ntemp,nsize,new_solid,rxgaos,rygaos,rzgaos,gndatoms,atoms_to_remove,ndatoms,&
                               rx,ry,rz)             
       
       use octahedre
       implicit none
       integer,intent(in)  :: ntemp,nsize,ndatoms,gndatoms
       integer, intent(out)  :: new_solid
       double precision,dimension(gndatoms),intent(in)     :: rxgaos,rygaos,rzgaos
       double precision,dimension(ndatoms),intent(in)     :: atoms_to_remove(ndatoms)
       double precision,dimension(nsize),intent(inout)     :: rx,ry,rz
       double precision,dimension(nsize)     :: rxtemp,rytemp,rztemp
       integer :: i,icnt,is,icard,i_test_card
       !
       !
       ! 
       icnt=0
       rxtemp(:)=99999999.d0
       rytemp(:)=99999999.d0
       rztemp(:)=99999999.d0
       write(*,*) ntemp, ndatoms
       do i=1,gndatoms
               icnt=icnt+1
             rxtemp(icnt)=rxgaos(i)
             rytemp(icnt)=rygaos(i)
             rztemp(icnt)=rzgaos(i)
                     write (77,'("atom Fe1   ",3f15.4)') rxgaos(icnt),rygaos(icnt),rzgaos(icnt)
      end do 
       
       
       
       i_test_card=0
       do i=1,ntemp
        icard=0
        do is=1,ndatoms
         if (i==nint(atoms_to_remove(is))) then
         i_test_card=i_test_card+1
         icard=1
         endif
      end do    
        if (icard==0) then
         icnt=icnt+1
         rxtemp(icnt)=rx(i)
         rytemp(icnt)=ry(i)
         rztemp(icnt)=rz(i)
        end if
      end do
        
        if (icnt.ne.nsize) then
         write(*,*) 'The vector size is not correct...icnt vs nsize:',icnt,nsize
         write(*,*) '<write_coordinates>'
       ! stop
       end if
        
        
         write(*,*) 'Final number of atoms.....: ',icnt
       new_solid=icnt
      
       rx(1:nsize)=rxtemp(1:nsize)!*volc
       ry(1:nsize)=rytemp(1:nsize)!*volc
       rz(1:nsize)=rztemp(1:nsize)!*volc

       return
       end


!#####################################################################
               
       subroutine write_out_octa (nall,ntemp,new_solid,gndatoms,volc,&
                       a0, rx, ry, rz, rxb, ryb, rzb, a, ndatoms,    &
                       atoms_to_remove, atoms_to_remove_type)
       
       implicit none
       
       integer  :: nall,ntemp,i,in,icnt,icnt2,gndatoms,new_solid,is
       double precision     :: rx(nall),ry(nall),rz(nall),a0,a02,CELLDM,ascl
       double precision, intent(in) :: rxb(nall),ryb(nall),rzb(nall)
       double precision     :: x,y,z, volc , delta       
       double precision,dimension(3,3) :: a
       integer, intent(in)            :: ndatoms
       double precision,intent(in)     :: atoms_to_remove(ndatoms)
       integer,intent(in)     :: atoms_to_remove_type(ndatoms)
       double precision :: xmin,ymin,zmin,xmax,ymax,zmax
       logical   :: lremove
  
       open(11,file='jsia.xyz')
       open(12,file='sia.xyz')
       open(19,file='usia.xyz')
       open(16,file='jusia.xyz')
       open(20,file='uvsia.xyz')
       open(21,file='uvbsia.xyz')
       open(23,file='uvbcsia.xyz')
       open(24,file='uvbcsia_sim.xyz')
       open(13,file='unit_cell',form='formatted') 
       open (17,file='geom_new.data',form='formatted',status='unknown')
       open (22,file='geom_org.data',form='formatted',status='unknown')
       open (18,file='formd_new_ndm_scl.str',form='formatted',& 
          status='unknown')
       
         CELLDM=1.0
       in=1
       
       a02=a0/2.d0
       ascl=MAXVAL(a*volc)

      open (15,file='formd_new.str',form='formatted',status='unknown')
        
      write(15,'(i7)') new_solid
        write(15,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(15,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(15,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


      write(17,'(i7)') new_solid
      write(17,*) 1.0
        write(17,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(17,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(17,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

      write(22,'(i7)') ntemp
      write(22,*) 1.0
        write(22,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(22,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(22,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


        write(18,'("1   1   1")') 
      write(18,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(18,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(18,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
      write(18,'(i8)') new_solid



10       format(I8,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
          write(13,*) 'number of atoms in the unit cell '
          write(13,*) new_solid 
          write(13,*) 'CELLDM= '
          write(13,*)  CELLDM
        write(13,*) 'position of the atoms in the unit cell'



       write(12,'(i7)') new_solid
       write(12,'(a)') ' '
       write(11,'(i7)') new_solid
       write(11,'(a)') ' '
       write(19,'(i7)') gndatoms
       write(19,'(a)') ' '
       write(16,'(i7)') gndatoms
       write(16,'(a)') ' '
       write(20,'(i7)') gndatoms+ndatoms
       write(20,'(a)') ' '

       write(21,'(i7)') gndatoms+ndatoms+ntemp
       write(21,'(a)') ' '



       xmin=1.d09 
       ymin=1.d09
       zmin=1.d09
                     
       xmax=-1.d09 
       ymax=-1.d09
       zmax=-1.d09

       do i=1,ntemp
             x = rxb(i)*volc 
             y = ryb(i)*volc 
             z = rzb(i)*volc 
        do is=1,ndatoms
         if (i==nint(atoms_to_remove(is))) then
         if (x .ge. xmax) xmax=x
         if (y .ge. ymax) ymax=y       
         if (z .ge. zmax) zmax=z
         
         if (x .le. xmin) xmin=x
         if (y .le. ymin) ymin=y       
         if (z .le. zmin) zmin=z
         
         endif
      end do           
        end do


       icnt2=0
       delta=+0.1
       xmin=xmin-delta
       ymin=ymin-delta
       zmin=zmin-delta
       xmax=xmax+delta
       ymax=ymax+delta
       zmax=zmax+delta

       do i=1,ntemp
         x = rxb(i)*volc 
         y = ryb(i)*volc 
         z = rzb(i)*volc 
         if (x.gt.xmin) then
       if (y.gt.ymin) then      
       if (z.gt.zmin) then
      
         if (x.lt.xmax) then
       if (y.lt.ymax) then      
       if (z.lt.zmax) then
           icnt2=icnt2+1
       end if
       end if
       end if
       end if
       end if
       end if    
       end do

       write(23,'(i7)') gndatoms+ndatoms+icnt2
       write(24,'(i7)') gndatoms+ndatoms+icnt2
       write(23,'(a)') ' '
       !gndatoms - atoms in interstitial position - label  as Fe
 
       icnt=0       
       do i=1,gndatoms
       icnt =icnt+1
         x = rx(i)*volc 
         y = ry(i)*volc 
         z = rz(i)*volc 
        write(12,'("Fe  ",3f18.13)') x,y,z
        write(11,'("Fe  ",3f18.13)') x,y,z
        write(19,'("Fe  ",3f18.13)') x,y,z
        write(16,'("Fe  ",3f18.13)') x,y,z
        write(20,'("Fe  ",3f18.13)') x,y,z
        write(21,'("Fe  ",3f18.13)') x,y,z
        write(23,'("Fe  ",3f18.13)') x,y,z
        write(24,'("Fe  ",3f18.13)') x,y,z
        write(13,10) icnt,x,y,z,in,in,in
        write(15,'(a3,3(f18.13))') 'Fe ', x,y,z
        write(17,'(a3,3(f18.13))') '   ', x,y,z
        write(18,'(a3,3(f18.13),"   1")') '   ', x/ascl,y/ascl,z/ascl
       end do




       do i=1,ntemp
             x = rxb(i)*volc 
             y = ryb(i)*volc 
             z = rzb(i)*volc 
            write(21,'("Cu  ",3f18.13)') x,y,z
        do is=1,ndatoms
         if (i==nint(atoms_to_remove(is))) then
         
           if (atoms_to_remove_type(is).eq.70) then
            write(20,'("He  ",3f18.13)') x,y,z
            write(21,'("He  ",3f18.13)') x,y,z
            write(23,'("He  ",3f18.13)') x,y,z
         end if 
           if (atoms_to_remove_type(is).eq.71) then
            write(20,'("H  ",3f18.13)') x,y,z
            write(21,'("H  ",3f18.13)') x,y,z
            write(23,'("H  ",3f18.13)') x,y,z
         end if 
         endif
      end do           
        end do

         write(*,*) xmin,ymin,zmin,xmax,ymax,zmax


       
       
       do i=gndatoms+1,new_solid
       write(22,'(a3,3(f18.13),"   1")')  '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
       icnt=icnt+1
       write(12,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc       
       write(11,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc       
       write(13,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
       write(15,'(a3,3(f18.13))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
       write(17,'(a3,3(f18.13))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
       write(18,'(a3,3(f18.13),"   1")') '   ', rx(i)*volc/ascl,        &
                                        ry(i)*volc/ascl, rz(i)*volc/ascl
       end do
       
      
       
       
       icnt=0
       do i=1,ntemp
        icnt =icnt+1
        x = rxb(i)*volc 
        y = ryb(i)*volc 
        z = rzb(i)*volc 
        if (x.gt.xmin) then
         if (y.gt.ymin) then      
          if (z.gt.zmin) then
           if (x.lt.xmax) then
            if (y.lt.ymax) then      
             if (z.lt.zmax) then
                write(23,'("Cu  ",3f18.13)') x,y,z
                lremove=.false.
                do is=1,ndatoms
                 if (i==nint(atoms_to_remove(is)))  lremove=.true.
              end do           
                 if (.not.lremove) write(24,'("Cu  ",3f18.13)') x,y,z
            
                

             end if
            end if
           end if
          end if
         end if
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
       


       return
       end
       


    subroutine write_formd_new_org (ntemp,nsize,rx,ry,rz,a,a0,volc)
    
    implicit none
    integer, intent(in)  :: ntemp,nsize
    real(8), intent(in)  :: rx(nsize),ry(nsize),rz(nsize)
    real(8), intent(in)  ::a(3,3),a0,volc
    real(8)   :: a02,ascl
    integer ::  i    
    
       
              a02=a0/2.d0
       ascl=MAXVAL(a*volc)
       open (21,file='formd_new_ndm_org.str',form='formatted',status='unknown')

        write(21,'("1   1   1")') 
      write(21,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(21,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(21,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
      write(21,'(i8)') ntemp

       do i=1,ntemp
       write(21,'("   ",3(f18.13),"   1")') rx(i)*volc/ascl,ry(i)*volc/ascl, rz(i)*volc/ascl
       end do
       
      close(21) 
    return
    end
