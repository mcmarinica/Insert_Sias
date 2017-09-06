    program put_vat
       use octahedre 
       implicit none
       
       real(8), parameter :: a0=2.8553,as110=2.05d0,as111=2.05d0,as100=2.01d0
       integer,dimension(:), allocatable :: n1,n2,n3,type_octa
       double precision,dimension(:),allocatable :: rx,ry,rz,rxb,ryb,rzb
       integer,dimension(:),allocatable :: ns1,ns2,ns3,nsite,indx,nsite_octa_vac_remove,nsite_octa_vac_remove_type,atoms_to_remove_type
       real(8),dimension(:),allocatable :: list_no_order,atoms_to_remove,offlattice_no_order,&
                                           rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos
       integer,dimension(:,:),allocatable :: nsite_m
       double precision     :: rsx(2),rsy(2),rsz(2)
       double precision     :: rgx(6),rgy(6),rgz(6)
       real(8), dimension(3,3) :: a
       integer :: tsia,nsia,nsia_a,ntemp,nall,nsize
       integer :: n1MAX,n2MAX,n3MAX,ndatoms,gndatoms,new_solid,nlac_remove,natoms_extra	
       real(8)  :: volc
       
       
       write(*,*) 'Choice the sia (1 - 110)(2 - 111)(3 - 100) (4-gao)(5 -vac)(6-octa):'
       read(*,*) tsia
       write(*,*) 'excelent choice master'
       
       if (tsia.gt.6)  stop
       
       call read_ntemp(ntemp)
       
         call read_nsia (nsia)
       
       if ( tsia.eq.4) then
         nsia_a=nsia*4
         nall  =ntemp+2*nsia
        else if (tsia.eq. 5) then
         nsia_a=-nsia
	 nall=ntemp-nsia
	else if (tsia.eq.6) then
	call init_octa
	
	nsia_a=nsia*10
	nall=ntemp + 2*nsia
       else
	 nsia_a=nsia
         nall  =ntemp+nsia
       end if	 
       
       print*, ntemp
       write(*,'("The number of lattice atoms .......:",i7)') ntemp
       write(*,'("The number of lines in no00.inp....:",i7)') nsia
       write(*,'("The number of no-lattice atoms ....:",i7)') nsia_a
       
       nsize=nall
       if (nall.lt.ntemp) nsize=ntemp
       allocate (n1(nsize),  n2(nsize),  n3(nsize),           &
                  rx(nsize),  ry(nsize), rz(nsize),           &
                  rxb(nsize),  ryb(nsize), rzb(nsize),           &
		ns1(nsia),ns2(nsia),ns3(nsia),nsite(nsia), &
		nsite_m(4,nsia),nsite_octa_vac_remove(10*nsia),nsite_octa_vac_remove_type(10*nsia),type_octa(nsia),indx(10*nsia), &
		list_no_order(10*nsia),offlattice_no_order(12*nsia), &
		temprxgaos(12*nsia),temprygaos(12*nsia),temprzgaos(12*nsia))

       !we read the input cell of the perfect bulk ...
       call read_input(ntemp,nsize,a,rx,ry,rz)
       
       
       if (nsia_a.gt.0) then
        call read_sia (tsia,as110,as111,as100,rsx,rsy,rsz,rgx,rgy,rgz)
       end if
       
       ! put the formd.str file in the integer n1,n2,n3...
       call write_n123 (ntemp,nsize,a0,n1,n2,n3,rx,ry,rz,  &
                             n1MAX,n2MAX,n3MAX)  
			     
       ! the the integer coordinates of the site where the inst must 
       ! be put
      if (tsia.ne.6) then      
       call read_config_sia (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
       type_octa(:)=1
      else if (tsia.eq.6) then
       write(*,*) 'IM HEEEEEEEEEEEEEEEERE'
       call read_config_sia_octa (nsia,type_octa,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
      end if	  
			     
       !
       call check_site (tsia,ntemp,nsize,nsia,nlac_remove,n1,n2,n3,ns1,ns2,ns3,&
                        nsite,nsite_m,nsite_octa_vac_remove,nsite_octa_vac_remove_type,type_octa)

       write(*,*) 'nlac_remove.........',nlac_remove
       volc=1.d0
       if (nsia_a.gt.0) then
       volc= (dble(nall)/dble(ntemp))**(1.d0/3.d0)
       end if
       print *, nall, ntemp, volc		         

       if (tsia.le.3) then
         call write_out   (nall,ntemp,nsia,volc,nsite,rx,ry, &
                       rz,rsx,rsy,rsz,a,a0)
		
       else if (tsia.eq.4) then
         call write_out_m (nall,ntemp,nsia,volc,a0,nsite_m,rx,ry, &
                       rz,rgx,rgy,rgz,a)

       else if (tsia.eq.5) then
       volc=1.d0
         call write_out_v (nall,ntemp,nsia,nsize,volc,nsite,rx,ry,rz,a)
       else if (tsia.eq.6) then	
         write(*,*) 'ntemp',ntemp
         write(*,*) 'nsia',nsia
         write(*,*) 'nsize',nsize
	 
	 call write_formd_new_org (ntemp,nsize,rx,ry,rz,a,a0,volc)
         rxb(:)=rx(:)
         ryb(:)=ry(:)
         rzb(:)=rz(:)
	 
	 call  lattice_compute_different_removeatoms (nsia, nlac_remove,nsite_octa_vac_remove,indx,list_no_order,ndatoms)
	 allocate (atoms_to_remove(ndatoms),atoms_to_remove_type(ndatoms))
	 call  lattice_compute_atoms_to_remove(atoms_to_remove,atoms_to_remove_type,ndatoms,nlac_remove,list_no_order,nsite_octa_vac_remove_type,indx)
	 call  offlattice_get_coordinates_and_order(a0,temprxgaos,temprygaos,temprzgaos,&
                                                    rx,    ry,    rz,    &
                     offlattice_no_order,indx,nsite,type_octa,nsia,nsize,gndatoms,natoms_extra)
		     
         allocate (rxgaos(gndatoms),rygaos(gndatoms),rzgaos(gndatoms))
	 call offlattice_compute_atoms_to_remove(gndatoms,rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos,&
	                       offlattice_no_order,indx,nsia,natoms_extra)

	 call write_coordinates (ntemp,nsize,new_solid,rxgaos,rygaos,rzgaos,gndatoms,atoms_to_remove,ndatoms,&
	                         rx,ry,rz)	       
		
        volc= (dble(new_solid)/dble(ntemp))**(1.d0/3.d0)
        write(*,*) 'The no of interstitials in defect ...:', ntemp-ndatoms+gndatoms 
	
	call  write_out_octa (nall,ntemp,new_solid,gndatoms,volc,a0,rx,  &
                       ry,rz,rxb,ryb,rzb,a,ndatoms,atoms_to_remove,atoms_to_remove_type)
       end if	 		

         call neighb_sia (nall,nsia,volc,nsite,rx,  &
                       ry,rz,a)

        end  program put_vat
	
!#####################################################################
!#####################################################################
