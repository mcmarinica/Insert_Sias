    program put_vct
       use math
       use octahedre
       use icosaedre
       use migration_sia 
       use the_111_smooth
       use the_100_smooth
       implicit none
       
       real(8), parameter :: a0=2.8553,as110=2.05d0,as111=1.60d0,as100=1.70d0
       integer,dimension(:), allocatable :: n1,n2,n3,type_octa
       real(8) ,dimension(:),allocatable :: rx,ry,rz,rxb,ryb,rzb
       integer,dimension(:),allocatable :: ns1,ns2,ns3,nsite,   &
                                           indx,nsite_octa_vac_remove,nsite_octa_vac_remove_type, &
                                            atoms_to_remove_type,isite_surface,isite_surfaceb
       real(8),dimension(:),allocatable :: list_no_order,atoms_to_remove,offlattice_no_order,&
                                           rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos
       integer,dimension(:,:),allocatable :: nsite_m
       double precision     :: rsx(2),rsy(2),rsz(2)
       double precision     :: rgx(6),rgy(6),rgz(6)
       real(8), dimension(3,3) :: a, amat, inv_amat,abase, inv_abase
       integer :: tsia
       integer :: nsia,nsia_a,ntemp,nall,nsize
       integer :: is_tot, is_totb, i,j,pdist, izcoordinate, n1MAX,n2MAX,n3MAX,ndatoms,gndatoms,new_solid,nlac_remove,natoms_extra 
       real(8)  :: volc
       real(8) , dimension(3,4) :: plane_corner
       real(8)                  :: surface
       real(8), dimension(3)    :: shift,temp, ctemp 
       real(8) :: dist_current, dist_ini
       integer ::  it_a,is,i_count
       !tsia - is integer type of sia
       call read_tsia (tsia)
       !ntemp number of atoms in reference structure (aka bulk)
       call read_ntemp(ntemp)
       !read nsia - number of defects the header of no00.inp or no00m.inp file
       call read_nsia (nsia,tsia)
       
       if ( tsia.eq.4) then
         nsia_a=nsia*4
         nall  =ntemp+2*nsia
        else if (tsia.eq.5) then
         nsia_a=-nsia
         nall=ntemp-nsia
        else if ((tsia.eq.6).or.(tsia==67)) then
         call init_octa
         nsia_a=nsia*10
         nall=ntemp + 2*nsia_a ! I do not understand. ilogical
        else if (tsia.eq.7) then
         call init_icos
         nsia_a=13*nsia   ! one in center and 6 dumbells i.e. 12 off_latice atoms 
         nall = ntemp+nsia_a! I do nor understand, ilogical
       else
         nsia_a=nsia
         nall  =ntemp+nsia
       end if
       
       print*, ntemp
       write(*,'("The number of lattice atoms .......:",i8)') ntemp
       write(*,'("The number of lines in no00.inp....:",i8)') nsia
       write(*,'("The number of no-lattice atoms ....:",i8)') nsia_a
       
       nsize=nall
       if (nall.lt.ntemp) nsize=ntemp
       allocate ( n1(nsize),  n2(nsize),  n3(nsize),                  &
                  rx(nsize),  ry(nsize), rz(nsize),                   &
                  rxb(nsize),  ryb(nsize), rzb(nsize),                &
                  ns1(nsia),ns2(nsia),ns3(nsia),nsite(nsia),          &
                  nsite_m(4,nsia),nsite_octa_vac_remove(10*nsia),     &
                  nsite_octa_vac_remove_type(10*nsia),type_octa(nsia),&
                  indx(10*nsia), list_no_order(10*nsia),              &
                  offlattice_no_order(12*nsia), temprxgaos(12*nsia),  &
                  temprygaos(12*nsia),temprzgaos(12*nsia),            &
                  isite_surfaceb(nsize), isite_surface(nsize))

       !we read the input cell of the perfect bulk ...return rx,ry, zy and the box at (through pbc module)
       call read_input(ntemp,nsize,a,rx,ry,rz)
       
       
       if (nsia_a.gt.0) then
        call read_sia (tsia,as110,as111,as100,rsx,rsy,rsz,rgx,rgy,rgz)
        if (tsia==21) then
         !
          call allocate_the_111 (nsia)
         !
        end if 

        if (tsia==31) then
         !
          call allocate_the_100 (nsia)
         !
        end if 
       end if
       
       ! put the formd.str file in the integer n1,n2,n3...
       call write_n123 (ntemp,nsize,a0,n1,n2,n3,rx,ry,rz,  &
                             n1MAX,n2MAX,n3MAX)  
     
       ! the the integer coordinates of the site where the inst must 
       ! be put
      if ((tsia.le.6).or.(tsia==21).or.(tsia==31).or.(tsia==67).or.(tsia==7)) then      
        !check inconsistencies 
        call read_config_sia (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
        if (tsia==7) call fix_ico_type (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
        !some driver for dc model
        call analyze_config_sia (tsia,nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)

        if (tsia==21) call fill_the_111 (a0,as111,rsx,rsy,rsz)
        if (tsia==31) call fill_the_100 (a0,as100,rsx,rsy,rsz)

        if (tsia==67) then
          write(*,*) n1MAX
            if ((n1MAX/=n2MAX).or.(mod(n1MAX-1,4)/=0)) then
              ! n1MA1X-1 because we have 
              write(*,'("n1MAX-1 and n2MAX-1 should be equal and multiple of 4...:", 2i4)') n1MAX-1,n2MAX-1
              write(*,'("stop")')
              stop
            end if
        end if 
        if ((tsia==6).or.(tsia==67)) then 
          !read the type_octa.useful for C15 case because, there are, two orientations and consequently two type of octa
          type_octa(:)=1
          call read_config_sia_octa (nsia,type_octa,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
        end if
      end if 


      if ((tsia==22).or.(tsia==221)) then
        call allocate_migration(nsia)
        if (tsia==221) then
          call allocate_the_111_migration (ntemp,nsia)
        end if 
        call read_config_sia_migration (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX) 
        if (tsia==221)then
          call fill_the_111 (a0,as111,rsx,rsy,rsz)
        end if 
      end if     
       !
      if (tsia==21) then
        call build_the_nsite_the_111(nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
      end if 
      if (tsia==31) then
         call build_the_nsite_the_100(nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
      end if 


      if (tsia==221) then
         call  build_the_nsite_the_111_migration (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
      end if 
      !check site and many more
      call check_site (tsia,ntemp,nsize,nsia,nlac_remove,n1,n2,n3,ns1,ns2,ns3,&
                        nsite,nsite_m,nsite_octa_vac_remove,nsite_octa_vac_remove_type,type_octa)
      if ((tsia==22).or.(tsia==221)) then 
         call check_site_migration (ntemp,nsize,nsia,n1,n2,n3,ns1,ns2,ns3)
           if (tsia==221) then
            call care_check_site_migration (ntemp,nsize,nsia,n1,n2,n3)
           end if 
       end if 

       write(*,*) 'nlac_remove (before suppression).........:',nlac_remove
       volc=1.d0
       if (nsia_a.gt.0) then
       volc= (dble(nall)/dble(ntemp))**(1.d0/3.d0)
       end if
       volc=1.d0
       print *, nall, ntemp, volc                   

       if ((tsia.le.3).or.(tsia==21).or.(tsia==31)) then
         rxb(:)=rx(:)
         ryb(:)=ry(:)
         rzb(:)=rz(:)
       
         if (tsia==21) then
            call distord_lattice_111 (ntemp,nsia,rx,ry,rz)
         end if 
         if (tsia==31) then
            call distord_lattice_100 (ntemp,nsia,rx,ry,rz)
         end if 
         call write_out   (nall,ntemp,nsia,volc,nsite,rxb,ryb,rzb,rx,ry, &
                       rz,rsx,rsy,rsz,a,a0)
          
       else if (tsia.eq.4) then
         call write_out_m (nall,ntemp,nsia,volc,a0,nsite_m,rx,ry, &
                       rz,rgx,rgy,rgz,a)

       else if (tsia.eq.5) then
         volc=1.d0
         call write_out_vv (nall,ntemp,nsia,nsize,volc,nsite,rx,ry,rz,a)
       else if ((tsia==7).or.(tsia.eq.6).or.(tsia==67)) then     
         write(*,*) 'ntemp',ntemp
         write(*,*) 'nsia',nsia
         write(*,*) 'nsize',nsize
      
         call write_formd_new_org (ntemp,nsize,rx,ry,rz,a,a0,volc)
         rxb(:)=rx(:)
         ryb(:)=ry(:)
         rzb(:)=rz(:)
         ! redifine the atoms that should be removed. In the case that are 
         ! atoms_to_be_removed double defined (e.g. in the case
         ! this subroutine will define ndatoms the atoms that will be removed  
         call  lattice_compute_different_removeatoms (nsia, nlac_remove,nsite_octa_vac_remove,indx,list_no_order,ndatoms)
         !the atoms that should be removed are packed into atoms_to_remove vector
         allocate (atoms_to_remove(ndatoms),atoms_to_remove_type(ndatoms))
         call  lattice_compute_atoms_to_remove(atoms_to_remove,atoms_to_remove_type,ndatoms, & 
              nlac_remove,list_no_order,nsite_octa_vac_remove_type,indx)
         if ((tsia.eq.6) .or. (tsia.eq.67)) then
           call  offlattice_get_coordinates_and_order(a0,temprxgaos,temprygaos,temprzgaos,&
                                                    rx,    ry,    rz,    &
                     offlattice_no_order,indx,nsite,type_octa,nsia,nsize,gndatoms,natoms_extra)
           allocate (rxgaos(gndatoms),rygaos(gndatoms),rzgaos(gndatoms))
           call offlattice_compute_atoms_to_remove(gndatoms,rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos,&
                            offlattice_no_order,indx,nsia,natoms_extra)



         else if (tsia.eq.7) then
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

         call write_coordinates (ntemp,nsize,new_solid,rxgaos,rygaos,rzgaos,gndatoms,atoms_to_remove,ndatoms,&
                              rx,ry,rz)
     
         write(*,*)'new_solid.................................:', new_solid                            
         volc= (dble(new_solid)/dble(ntemp))**(1.d0/3.d0)
         volc=1.d0
         write(*,*)'The no of defects SIAs ...................:', -ndatoms+gndatoms 
     
        call  write_out_octa (nall,ntemp,new_solid,gndatoms,volc,a0,rx,  &
                       ry,rz,rxb,ryb,rzb,a,ndatoms,atoms_to_remove,      &
                       atoms_to_remove_type, nsite, nsia)
       
            !pdist should be 2*k number is we want a complete cubic cell
            !pdist=(n1MAX-1)/2
            !pdist=48  !10 nm with 48-ABAB 49-ABA ! with cube.x = 89a0
            !pdist=18  !3  nm ! with cube.x = 51a0 or 49a0 or 53a0
            pdist=25    !5 nm !    24-ABAB  25-ABA   with cube.x = 89a0
            !pdist=5
            !if pdist is 2*k+1 the cell is incomplete ...
      
            if  (tsia==67) then !Kazuto style
             !
             if ((pdist>(n1MAX-1)/2).or.(pdist<=0))  then
               write (*,'(" p should be < (n1MAX-1)/2; and > 0...:",2i7)') pdist,(n1MAX-1)/2
               write (*,*) 'stop in put_vct'
               stop
             end if 
             !   
             if (mod(pdist,2)/=0) then
               write (*,'(" p should be 2*k, k is an integer.....: ",2i7)') pdist,(n1MAX-1)/2
               !write (*,*) 'stop in put_vct'
               !stop
             end if
               
             call define_center (pdist,n1MAX,surface,plane_corner)
             izcoordinate=nint(dble(n1MAX-1)/2.d0*dsqrt(2.d0))
             if (mod(izcoordinate,2)/=0) then
               izcoordinate=izcoordinate+1
             end if 
             izcoordinate = izcoordinate+1
             if (izcoordinate.gt.(n3MAX-1)) then
              write(*,*) 'At least the n3MAX should be increased ....', n3MAX
               stop
             end if 

              
            call solid_and_defect_cut (a0, surface, new_solid, gndatoms, nsize, rx,ry,rz, isite_surface, &
                                       is_tot,izcoordinate, plane_corner)  
            call solid_and_defect_cut (a0, surface, ntemp, gndatoms, nsize, rxb,ryb,rzb, isite_surfaceb, &
                                       is_totb,izcoordinate, plane_corner)  
             open (unit=44, file='k_bulk.gin',status='unknown')
             open (unit=45, file='k_dfct.gin',status='unknown')
             open (unit=46, file='k_dfct_vide.gin',status='unknown')
             write(46,*) '1 1 1'
             write(45,*) '1 1 1'
             write(44,*) '1 1 1'
             
             shift(:)=(/ 1.d0,1.d0,0.d0 /)
             abase(:,1) = (a0/2.d0)*(plane_corner(:,2)-plane_corner(:,1) + shift(:))
             shift(:) = (/ -1.d0,1.d0,0.d0 /)
             abase(:,2) = (a0/2.d0)*(plane_corner(:,4)-plane_corner(:,1) + shift(:))
             abase(:,3) =  (a0/2.d0)*(/ 0.d0, 0.d0, dble(izcoordinate)+1.d0 /)

             abase(:,1) = (/ dsqrt(DOT_PRODUCT(abase(:,1),abase(:,1))), 0.d0, 0.d0 /)
             abase(:,2) = (/ 0.d0, dsqrt(DOT_PRODUCT(abase(:,2),abase(:,2))),       0.d0 /)
             abase(:,3) = (/ 0.d0, 0.d0, dsqrt(DOT_PRODUCT(abase(:,3),abase(:,3)))  /)



 
             call MatInv(abase,inv_abase)
             write(45,'(3f15.8)') abase(:,1)
             write(45,'(3f15.8)') abase(:,2)  
             write(45,'(3f15.8)') abase(:,3) 
             write(45,'(i6)') is_tot

             write(44,'(3f15.8)') abase(:,1)
             write(44,'(3f15.8)') abase(:,2)  
             write(44,'(3f15.8)') abase(:,3) 
             write(44,'(i6)') is_totb


             write(46,'(3f15.8)') abase(:,1)
             write(46,'(3f15.8)') abase(:,2)  
             write(46,'(3f15.8)') abase(:,3) 
             write(46,'(i6)') gndatoms+nsia



             amat(1,:) = (/  1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0),  0.d0 /) 
             amat(2,:) = (/ -1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0),  0.d0 /) 
             amat(3,:) = (/             0.d0,            0.d0,  1.d0 /) 

             call MatInv(amat,inv_amat)

             shift(:) = plane_corner(:,1)*(a0/2.d0)

             write(*,'("TTTTTTTTTTTT", 3f12.6)') amat(1,1)*abase(1,1)/(a0/2.d0), amat(2,1)*abase(2,2)/(a0/2.d0)
             
             !The suggested C15 center:
              do j=1,3
               temp(j)=shift(j) + amat(1,j)*abase(1,1)/2.d0+ amat(2,j)*abase(2,2)/2.d0+amat(3,j)*abase(3,3)/2.d0
              end do
              write(*,'("The suggested center for C15.....:",3i6)') int(temp(:)/(a0/2.d0))

              
             
             do i=1,gndatoms
                  temp(:) = (/ rx(i), ry(i) , rz(i) /)
                  !transfer into the new  oXYZ 
                  do j=1,3 
                   ctemp(j) =DOT_PRODUCT(inv_amat(:,j),temp(:)-shift(:))
                  end do

                  if (isite_surface(i)==1) then
                  !    write(55,'(i6,a3,3f10.4,2f7.2)')i,' Fe',ctemp(1:3),1.d0,0.5d0
                    write(45,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
                    write(46,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
                  end if 
            end do


            do is=1,nsia
                  i=nsite(is)
                  temp(:) = (/ rxb(i), ryb(i) , rzb(i) /)
                  !transfer into the new  oXYZ 
                  do j=1,3 
                   ctemp(j) =DOT_PRODUCT(inv_amat(:,j),temp(:)-shift(:))
                  end do

                  !    write(55,'(i6,a3,3f10.4,2f7.2)')i,' Fe',ctemp(1:3),1.d0,0.5d0
                    write(45,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
                    write(46,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
            end do




      
             i_count=0
             do i=gndatoms+1,new_solid

               dist_ini=99999999.9
                !
                do is=1,nsia
                  it_a=nsite(is)
                  dist_current=dsqrt((rx(i)-rxb(it_a))**2+(ry(i)-ryb(it_a))**2+(rz(i)-rzb(it_a))**2)
                  if (dist_current < dist_ini) dist_ini=dist_current
                end do
                !

               if (dist_ini  > 0.01) then
                  i_count=i_count+1
                  temp(:) = (/ rx(i), ry(i) , rz(i) /)
                  !transfer into the new  oXYZ 
                  do j=1,3 
                   ctemp(j) =DOT_PRODUCT(inv_amat(:,j),temp(:)-shift(:))
                  end do

                  if (isite_surface(i)==1) then
                  !    write(55,'(i6,a3,3f10.4,2f7.2)')i,' Fe',ctemp(1:3),1.d0,0.5d0
                    write(45,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
                  end if 
                end if 

            end do

!debug            write(*,*) gndatoms, nsia, i_count
!debug            stop




             do i=1,ntemp
                temp(:) = (/ rxb(i), ryb(i) , rzb(i) /)
                !transfer into the new reper 
                do j=1,3 
                  ctemp(j) =DOT_PRODUCT(inv_amat(:,j),temp(:)-shift(:))
                end do
                !
                if (isite_surfaceb(i)==1) then
                    write(44,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
                end if 
            end do
            close (44)
            close (45) 

            end if    ! tsia==67  
             
             



      else if ((tsia == 22).or.(tsia == 221)) then
      rxb(:)=rx(:)
      ryb(:)=ry(:)
      rzb(:)=rz(:)

!    I should change here .... but there are two cases deb and fin ....
!        call write_out   (nall,ntemp,nsia,volc,nsite,rx,ry, &
!                       rz,rsx,rsy,rsz,a,a0)
       if (tsia==221) then
         call deb_distord_lattice_111_migration (ntemp,nsia,rxb,ryb,rzb,rx,ry,rz)
       end if 
       call deb_fill_the_mig_objects (ntemp,rx,ry,rz,rsx,rsy,rsz)
       call deb_write_out_migration   (nall,ntemp,volc,rx,ry, &
                       rz,rsx,rsy,rsz,a)

       call write_out   (nall,ntemp,nsia,volc,nsite,rxb,ryb,rzb,rx,ry, &
                       rz,rsx,rsy,rsz,a,a0)



       if (tsia==221) then
           call fin_distord_lattice_111_migration (ntemp,nsia,rxb,ryb,rzb,rx,ry,rz)
       end if 
       call fin_fill_the_mig_objects (ntemp,rx,ry,rz,rsx,rsy,rsz)
       call fin_write_out_migration   (nall,ntemp,volc,rx,ry, &
                       rz,rsx,rsy,rsz,a)
 

       end if  !end of big else if          

    !call neighb_sia (nall,nsia,volc,nsite,rx,  &
    !                   ry,rz,a)

        end  program put_vct
!#####################################################################
!#####################################################################
