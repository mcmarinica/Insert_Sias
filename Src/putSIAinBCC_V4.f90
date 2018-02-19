

       implicit double precision (a-h,o-z)
       parameter (a0=2.8553,as110=2.05d0,as111=2.05d0,as100=2.05d0)
       integer,dimension(:), allocatable :: n1(:),n2(:),n3(:)
       double precision,dimension(:),allocatable :: rx(:),ry(:),rz(:)
       integer,dimension(:),allocatable :: ns1(:),ns2(:),ns3(:),nsite(:)
       integer,dimension(:,:),allocatable :: nsite_m
       double precision     :: rsx(2),rsy(2),rsz(2)
       double precision     :: rgx(6),rgy(6),rgz(6)
       double precision,dimension(:,:) :: a(3,3)
       integer :: tsia,nsia
       
       
       write(*,*) 'Choice the sia (1 - 110)(2 - 111)(3 - 100) (4-gao)(5 -vac):'
       read(*,*) tsia
       write(*,*) 'excelent choice master'
       
       if (tsia.gt.5)  stop
       
       call read_ntemp(ntemp)
       
       call read_nsia (nsia)
       
       if ( tsia.eq.4) then
         nsia_a=nsia*4
         nall  =ntemp+2*nsia
        else if (tsia.eq. 5) then
         nsia_a=-nsia
	 nall=ntemp-nsia
       else
	 nsia_a=nsia
         nall  =ntemp+nsia
       end if	 
       
       print*, ntemp
       print*, nsia, nsia_a

       
       nsize=nall
       if (nall.lt.ntemp) nsize=ntemp
       allocate (n1(nsize),  n2(nsize),  n3(nsize),           &
                  rx(nsize),  ry(nsize), rz(nsize),           &
		ns1(nsia),ns2(nsia),ns3(nsia),nsite(nsia), &
		nsite_m(4,nsia))

       !we read the input cell ...
       call read_input(nall,ntemp,nsize,a,rx,ry,rz)
       
       
       if (nsia_a.gt.0) then
        call read_sia (tsia,as110,as111,as100,rsx,rsy,rsz,rgx,rgy,rgz)
       end if
       
       ! put the formd.str file in the integer n1,n2,n3...
       call write_n123 (nall,ntemp,nsize,a0,a,n1,n2,n3,rx,ry,rz,  &
                             n1MAX,n2MAX,n3MAX)  
			     
       ! the the integer coordinates of the site where the inst must 
       ! be put
       call read_config_sia (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
			     
       !
       call check_site (tsia,nall,ntemp,nsize,nsia,n1,n2,n3,ns1,ns2,ns3,&
                        nsite,nsite_m)

       
       volc=1.d0
       if (nsia_a.gt.0) then
       volc= (dble(nall)/dble(ntemp))**(1.d0/3.d0)
       end if
       print *, nall, ntemp, volc		         

       if (tsia.le.3) then
         call write_out (nall,ntemp,nsia,volc,n1,n2,n3,nsite,rx,ry, &
                       rz,rsx,rsy,rsz,a)
		
       else if (tsia.eq.4) then
         call write_out_m (nall,ntemp,nsia,volc,a0,n1,n2,n3,nsite_m,rx,ry, &
                       rz,rgx,rgy,rgz,a)

       else if (tsia.eq.5) then
         call write_out_v (nall,ntemp,nsia,nsize,volc,n1,n2,n3,nsite,rx,ry, &
                       rz,rsx,rsy,rsz,a)
		
       end if	 		

         call neighb_sia (nall,ntemp,nsia,volc,n1,n2,n3,nsite,rx,  &
                       ry,rz,rsx,rsy,rsz,a)

        end 
	
!#####################################################################
!#####################################################################
!#####################################################################
	 	
	
       subroutine read_ntemp (ntemp)
       integer :: ntemp
        
	open (10,file='formd.str',form='formatted',status='unknown')
        read (10,'(i8)') ntemp
       
       return
       end
!#####################################################################
	
       subroutine read_nsia (nsia)
       integer :: nsia
        
	open (11,file='no00.inp',form='formatted',status='unknown')
        read (11,*) nsia
       
       return
       end
!#####################################################################
       	  
       subroutine read_input(nall,ntemp,nsize,a,rx,ry,rz)      
       integer              :: nall, ntemp
       double precision     :: rx(nsize),ry(nsize),rz(nsize)
       double precision,dimension(:,:) :: a(3,3)
       character*3 :: sym
       

        read(10,'(3f18.13)')a(1,1), a(1,2), a(1,3)
        read(10,'(3f18.13)')a(2,1), a(2,2), a(2,3)
        read(10,'(3f18.13)')a(3,1), a(3,2), a(3,3)
	
	do i=1,ntemp
		read(10,'(a3,3(f18.13))') sym, rx(i),ry(i),rz(i)
        end do 		
     
		
       
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
	 else if (tsia==2) then
          as2 = as111/2.d0
          rsx(1)= - cos01 * cos45 *  as2
	  rsy(1)= - cos01 * cos45 *  as2
	  rsz(1)= - cos01 * cos45 *  as2
	  
	  rsx(2)= cos01 * cos45 *  as2
	  rsy(2)= cos01 * cos45 *  as2
	  rsz(2)= cos01 * cos45 *  as2
	 else if (tsia==3) then
          rsx(1)= 0.d0
	  rsy(1)= 0.d0
	  rsz(1)= -as100/2.d0
	  
	  rsx(2)= 0.d0
	  rsy(2)= 0.d0
	  rsz(2)= as100/2.d0
	  else if (tsia==4) then
	  
	  open (23,file='gao_unit.in')
	   
	  do i=1,6 
	   read(23,*) rgx(i),rgy(i),rgz(i)
	  end do
	  close(23)
	 else
	   print*, 'Unknown type of SIA: Check read_sia' 
	   stop
	 end if    
		
        
       return
       end 
!#####################################################################
       	  
       
       
       subroutine write_n123(nall,ntemp,nsize,a0,a,n1,n2,n3,rx,ry,rz,  &
                             n1MAX,n2MAX,n3MAX)      
       implicit double precision (a-h,o-z)
       integer  :: n1(nsize),n2(nsize),n3(nsize),i,n1MAX,n2MAX,n3MAX
       double precision     :: rx(nsize),ry(nsize),rz(nsize)
       double precision     :: a(3,3)
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
	  
	  print*, n1MAX,n2MAX,n3MAX
       
       return
       end 
!#####################################################################
       	  
       
       subroutine read_config_sia (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
       
       integer  :: nsia, n1MAX,n2MAX,n3MAX  
       integer  :: ns1(nsia),ns2(nsia),ns3(nsia)
            
	    
	    do i=1,nsia
	     read (11,*) ns1(i),ns2(i),ns3(i)
	     !write(*,*) 'ns11',ns1(i),ns2(i),ns3(i)
	      if ( (ns1(i) > n1MAX) .or.(ns2(i) > n2MAX) .or. &
	        (ns3(i) > n3MAX) ) then
		print*, 'Big problem in read_config_sia'
		print*, 'The cell is not big enough. Increase the size!'
		stop
	      end if	    
	    end do 
	          
       return
       end 
!#####################################################################
       	  
       	  
       subroutine check_site (tsia,nall,ntemp,nsize,nsia,n1,n2,n3,ns1,ns2,ns3, &
                              nsite,nsite_m)		         
       
       implicit none
       integer :: nall,ntemp,itst,ii,i,tsia,is,nsia,nsize
       integer, dimension(:)  :: n1(nsize),n2(nsize),n3(nsize), &
                        nsite(nsia),ns1(nsia),ns2(nsia),ns3(nsia)
       integer, dimension(:,:) :: nsite_m(4,nsia),           &
                        ns1_m(4,nsia),ns2_m(4,nsia),ns3_m(4,nsia)		
				



	itst=0
	do is=1,nsia
	  write(*,*) 'ns1',ns1(1:3)
	   do i=1,ntemp
	       if (n1(i)==ns1(is)) then
	        if (n2(i)==ns2(is)) then
	         if (n3(i)==ns3(is)) then
		   itst =itst +1
		   nsite(is)= i
		 end if 
		end if
	       end if
	    end do
	  end do     			

        if (itst.ne.nsia) then
	  print*, 'there is a big problem in input'
	  print*, 'Check check_site'
	  print*, itst, nsia
	  stop  	   	
	end if	


       if (tsia.eq.4) then	
        
	

	do ii=1,nsia
	   nsite_m(1,ii)=nsite(ii)

	   ns1_m(1,ii) = ns1(ii)	
	   ns1_m(2,ii) = ns1(ii)+1	
           ns1_m(3,ii) = ns1(ii)-1	
           ns1_m(4,ii) = ns1(ii)+1	
	
           ns2_m(1,ii) = ns2(ii)	
	   ns2_m(2,ii) = ns2(ii)-1
           ns2_m(3,ii) = ns2(ii)-1
           ns2_m(4,ii) = ns2(ii)+1

           ns3_m(1,ii) = ns3(ii)	
           ns3_m(2,ii) = ns3(ii)-1
           ns3_m(3,ii) = ns3(ii)+1
           ns3_m(4,ii) = ns3(ii)+1
	end do
	
	
	itst=0
	do is=1,nsia
	   do i=1,ntemp
	    do ii=2,4
	       if (n1(i)==ns1_m(ii,is)) then
	        if (n2(i)==ns2_m(ii,is)) then
	         if (n3(i)==ns3_m(ii,is)) then
		   itst =itst +1
		   nsite_m(ii,is)= i
		 end if 
		end if
	       end if
	     end do  
	    end do
	  end do  
	     			
        
	if (itst.ne.3*nsia) then
	  print*, 'there is a big problem in input'
	  print*, 'Check check_site 3*'
	  print*, itst, nsia
	  stop  	   	
  	end if	
	
	end if  
				
	return
	end			       
!#####################################################################
       	  
       subroutine write_out (nall,ntemp,nsia,volc,n1,n2,n3,nsite,rx,  &
                       ry,rz,rsx,rsy,rsz,a)
       
       use math 
       use pbc
       integer  :: nall,ntemp,nsia,n,i,icard
       double precision     :: rx(nall),ry(nall),rz(nall)
       integer  :: n1(nall),n2(nall),n3(nall),nsite(nsia)       
       double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc        
       double precision,dimension(:,:) :: a(3,3)
       real(kind(1.d90)), dimension(3) :: xp,xc
       open(12,file='sia.xyz')
       open(19,file='usia.xyz')
       open(13,file='unit_cell',form='formatted') 
       open (17,file='geom_new.data',form='formatted',status='unknown')
       open (22,file='geom_org.data',form='formatted',status='unknown')
       open (18,file='formd_new_ndm_scl.str',form='formatted',status='unknown')
       open (21,file='formd_new_ndm_org.str',form='formatted',status='unknown')
       
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


10       format(I5,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
          write(13,*) 'number of atoms in the unit cell '
          write(13,*) nall 
          write(13,*) 'CELLDM= '
          write(13,*)  CELLDM
	  write(13,*) 'position of the atoms in the unit cell'



       write(12,'(i8)') nall
       write(12,'(a)') ' '
       write(19,'(i8)') 2*nsia
       write(19,'(a)') ' '
       icnt=0
       do i=1,nsia
        do n=1,2
	 icnt =icnt+1
         x = rx(nsite(i))*volc + rsx(n)*volc
         y = ry(nsite(i))*volc + rsy(n)*volc
         z = rz(nsite(i))*volc + rsz(n)*volc
       write(12,'("Fe  ",3f18.13)') x,y,z
       write(19,'("Fe  ",3f18.13)') x,y,z
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
	   if (i==nsite(is)) then
	   icard=1
	   endif
	end do    
       if (icard==0) then
       icnt=icnt+1
       write(12,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc       
       write(13,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
       write(15,'(a3,3(f18.13))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
       write(17,'(a3,3(f18.13))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
        xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
        xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
        write(18,'(a3,3(f18.13),"   1")') '   ' ,xc(1:3)
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
       close(13)
       close(14)
       


       return
       end
!#####################################################################
       
!#####################################################################
       	  
       subroutine write_out_m (nall,ntemp,nsia,volc,a0,n1,n2,n3,nsite_m,rx,  &
                       ry,rz,rsx,rsy,rsz,a)
       use math
       use pbc
       integer  :: nall,ntemp,nsia,n,i,icard
       double precision     :: rx(nall),ry(nall),rz(nall),a0,a02
       integer  :: n1(nall),n2(nall),n3(nall),nsite(nsia),nsite_m(4,nsia)       
       double precision     :: rsx(6),rsy(6),rsz(6),x,y,z, volc        
       double precision,dimension(:,:) :: a(3,3)
       real(kind(1.d0)), dimension(3) :: xc,xp 
       open(12,file='sia.xyz')
       open(19,file='usia.xyz')
       open(13,file='unit_cell',form='formatted') 
       open (17,file='geom_new.data',form='formatted',status='unknown')
       open (22,file='geom_org.data',form='formatted',status='unknown')
       open (18,file='formd_new_ndm_scl.str',form='formatted',status='unknown')
       open (21,file='formd_new_ndm_org.str',form='formatted',status='unknown')
       
         CELLDM=1.0
	 in=1
	 
	 a02=a0/2.d0

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


10       format(I5,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
          write(13,*) 'number of atoms in the unit cell '
          write(13,*) nall 
          write(13,*) 'CELLDM= '
          write(13,*)  CELLDM
	  write(13,*) 'position of the atoms in the unit cell'



       write(12,'(i8)') nall
       write(12,'(a)') ' '
       write(19,'(i8)') 6*nsia
       write(19,'(a)') ' '
       icnt=0
       do i=1,nsia
        do n=1,6
	 icnt =icnt+1
         x = rx(nsite_m(1,i))*volc + rsx(n)*volc*a0
         y = ry(nsite_m(1,i))*volc + rsy(n)*volc*a0
         z = rz(nsite_m(1,i))*volc + rsz(n)*volc*a0
       write(12,'("Fe  ",3f18.13)') x,y,z
       write(19,'("Fe  ",3f18.13)') x,y,z
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
       write(12,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc       
       write(13,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
       write(15,'(a3,3(f18.13))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
       write(17,'(a3,3(f18.13))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
        xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
        xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
        write(18,'(a3,3(f18.13),"   1")') '   ' ,xc(1:3)
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
       close(13)
       close(14)
       


       return
       end
!#####################################################################
       	  
       subroutine write_out_v (nall,ntemp,nsia,nsize,volc,n1,n2,n3,nsite,rx,  &
                       ry,rz,rsx,rsy,rsz,a)
       use math
       use pbc       
       integer  :: nall,ntemp,nsia,n,i,icard,nsize
       double precision     :: rx(nsize),ry(nsize),rz(nsize)
       double precision     :: rsx(2),rsy(2),rsz(2)
       integer  :: n1(nsize),n2(nsize),n3(nsize),nsite(nsia)       
       double precision     :: x,y,z, volc        
       double precision,dimension(:,:) :: a(3,3)
       real(kind(1.d0)), dimension(3) :: xc,xp

       open(12,file='sia.xyz')
       open(19,file='usia.xyz')
       open(13,file='unit_cell',form='formatted') 
       open (17,file='geom_new.data',form='formatted',status='unknown')
       open (22,file='geom_org.data',form='formatted',status='unknown')
       open (18,file='formd_new_ndm_scl.str',form='formatted',status='unknown')
       open (21,file='formd_new_ndm_org.str',form='formatted',status='unknown')
       
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


10       format(I5,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
          write(13,*) 'number of atoms in the unit cell '
          write(13,*) nall 
          write(13,*) 'CELLDM= '
          write(13,*)  CELLDM
	  write(13,*) 'position of the atoms in the unit cell'



       write(12,'(i8)') ntemp
       write(12,'(a)') ' '
       write(19,'(i8)') nsia
       write(19,'(a)') ' '
       icnt=0
       do i=1,nsia
         x = rx(nsite(i))*volc 
         y = ry(nsite(i))*volc 
         z = rz(nsite(i))*volc 
        write(12,'("Fe  ",3f18.13)') x,y,z
        write(19,'("Fe  ",3f18.13)') x,y,z
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
       write(12,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc       
       write(13,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
       write(15,'(a3,3(f18.13))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
       write(17,'(a3,3(f18.13))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
        xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
        xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
        write(18,'(a3,3(f18.13),"   1")') '   ' ,xc(1:3)
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
       close(13)
       close(14)
       


       return
       end
!#####################################################################
!#####################################################################

       subroutine neighb_sia (nall,ntemp,nsia,volc,n1,n2,n3,nsite,rx,  &
                       ry,rz,rsx,rsy,rsz,a)
       
       parameter (NVMAX=6)
       integer  :: nall,ntemp,nsia,n,i,icard,nsia2
       double precision     :: rx(nall),ry(nall),rz(nall)
       double precision     :: rs(3,nsia),ds(nsia*nsia),dsn(nsia*nsia)
       integer  :: n1(nall),n2(nall),n3(nall),nsite(nsia),indx(nsia*nsia)       
       double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc        
       double precision,dimension(:,:) :: a(3,3)
        integer          :: nvv(NVMAX), nvvf(NVMAX)
	double precision :: dv(NVMAX),dvf(NVMAX),bcc(NVMAX-1)        
      
         CELLDM=1.0
	 

        bcc(1)=1.d0
	bcc(2)=2.d0/dsqrt(3.d0)
	bcc(3)=dsqrt(2.d0)*2.d0/dsqrt(3.d0)
	bcc(4)=dsqrt(11.d0)/2.d0*bcc(2)
	bcc(5)=2.d0
	
	a(:,:)=a(:,:)*volc
	
        do i=1,nsia
	 rs(1,i)=rx(nsite(i))*volc
	 rs(2,i)=ry(nsite(i))*volc
	 rs(3,i)=rz(nsite(i))*volc
	end do

         
       call d_pbc(a,rs,ds,nsia)
       nsia2=nsia*nsia
       call indexx(nsia2,ds,indx) 
       do i=1,nsia2
        j=indx(i)
        dsn(i)=ds(j)
       end do 
       
       call countn(nsia2,NVMAX,dsn,nvv,dv)
       
      
      nvvf(:)=0
      nvvf(1)=nvv(1)
      do i=1,NVMAX-1
       bn=bcc(i)*dv(2)
       write(*,*) bn
       do j = 2,NVMAX
         if (dabs(dv(j)-bn).lt.0.001) then
	  nvvf(i+1)=nvv(j)
	  dvf(i+1)=dv(j)
	 end if
       end do
      end do      
       
 	do i=1,NVMAX
 	write(*,'(2i7,f12.5)') i,nvv(i),dv(i)
 	end do

       write(*,'("=========================================")')
       write(*,'("      n1","      n2","      n3","      n4","      n5","     ntot")')
       write(*,'("=========================================")')
       write(*,'(6i8)')   (nvvf(i)/2,i=2,6),nvvf(1) 
       write(*,'(6f8.3)') (dvf(i),i=2,6   ),dvf(1)
       return
       end
!#####################################################################
       
       
      subroutine countn(NRMAX,NVMAX,a,nvv,dv)
! input the vector a(NRMAX)
! the output: the different values in dv(NVMAX)
!             how many values are for dv(i) recorded in nvv(i)     
      implicit none
      integer :: NRMAX,NVMAX
      double precision  :: at,ac,  a(NRMAX),dv(NVMAX)
      integer ::  nvv(NVMAX),ij,i,iv
      	     
      at=a(1)
      ij=1
      iv=1
        do i=2,NRMAX
	  ac=dabs(a(i)-at)
	  if (ac.lt.0.01) then
	   ij=ij+1
	  else
	   iv=iv+1   
            if (iv.gt.NVMAX) then
	     go to 113
	    end if 
	   ij=1
	  end if
	   nvv(iv)=ij
	   dv(iv)=a(i)

  	write(*,'(2f8.2,2i8,f8.3)') a(i),ac,ij,iv,dv(iv)
	  at=a(i)
	 end do
	 
113     continue 
       return
       end     


!*********************************************************************!
        SUBROUTINE  d_pbc (unit_vect,tau,ds,nat_up)
        implicit none
        integer  nat_up,i,j,ni
	double precision x,y,z,a_side,b_side,c_side,a_side2,&
	                 b_side2,c_side2,a,b,c
        double precision tau(3,nat_up),     &
                          unit_vect(3,3),  &
			  dc(3,3),dci(3,3), &
			  ds( nat_up*nat_up )
!
!       determination du nombre de voisins jusqu'a Rmax
!	
	
!..............pbc setup

	a_side=sqrt(DOT_PRODUCT(unit_vect(:,1),unit_vect(:,1)))
	b_side=sqrt(DOT_PRODUCT(unit_vect(:,2),unit_vect(:,2)))
	c_side=sqrt(DOT_PRODUCT(unit_vect(:,3),unit_vect(:,3)))
	
	
	a_side2=a_side/2.d0
	b_side2=b_side/2.d0
	c_side2=c_side/2.d0
	
	dc(:,1)=unit_vect(:,1)/a_side
	dc(:,2)=unit_vect(:,2)/b_side
	dc(:,3)=unit_vect(:,3)/c_side
       
        call inverse(dc,dci)   

        ni=0
	Do  i=1,nat_up
         Do j=1,nat_up
            ni=ni+1
	    x = tau(1,i)-tau(1,j)          
            y = tau(2,i)-tau(2,j)          
            z = tau(3,i)-tau(3,j)   
             

            call fromto(x,y,z,a,b,c,1,dc,dci)

               if (a.lt.-a_side2    ) a=a+a_side
               if (a.gt. a_side2    ) a=a-a_side

               if (b.lt.-b_side2    ) b=b+b_side
               if (b.gt. b_side2    ) b=b-b_side

               if (c.lt.-c_side2    ) c=c+c_side
               if (c.gt. c_side2    ) c=c-c_side

            call fromto(x,y,z,a,b,c,-1,dc,dci)
    
           ds(ni)=dsqrt(x**2+y**2+z**2)


        end do
       end do
	   

       RETURN
       END
!*********************************************************************!


      subroutine fromto(x,y,z,a,b,c,imode,dc,dci)
!**********************************************************************
!*								      *
!* NAME: "FROMTO"						      *
!*								      *
!* USAGE: Transformation between crystallographic and Cartesian       *
!*	  coordinate frames.					      *
!*								      *
!* METHOD:			-1				      *
!*		[a]    [ax bx cx]  [x]       [x]		      *
!*		[b] =  [ay by cy]  [y] = DCI [y]		      *
!*		[c]    [az bz cz]  [z]       [z]		      *
!*								      *
!*		[x]    [ax bx cx]  [a]       [a]		      *
!*		[y] =  [ay by cy]  [b] = DC  [b]		      *
!*		[z]    [az bz cz]  [c]       [c]		      *
!*								      *
!* where matrix DC is a matrix of direction cosines of  	      *
!* crystallographic axes.					      *
!*								      *
!* INPUT/OUTPUT: 'x, y, z' - Cartesian coordinates;		      *
!*		 'a, b, c' - crystallographic coordinates;	      *
!*               'imode'   - direction of transformation:             *
!*                imode= 1  Cartesian --> crystallographic;           *
!*                imode=-1  crystallographic --> Cartesian.           *
!*                                                                    *                                                                    *
!* NOTE: Should be used together with (and after) MDCELL and MDCELL2. *
!**********************************************************************
      implicit double precision (a-h,o-z)
      double precision :: dc(3,3),dci(3,3)

      if (imode.eq.1) then

         a=dci(1,1)*x + dci(1,2)*y + dci(1,3)*z
         b=dci(2,1)*x + dci(2,2)*y + dci(2,3)*z
         c=dci(3,1)*x + dci(3,2)*y + dci(3,3)*z

      else if (imode.eq.-1) then

         x=dc(1,1)*a + dc(1,2)*b + dc(1,3)*c
         y=dc(2,1)*a + dc(2,2)*b + dc(2,3)*c
         z=dc(3,1)*a + dc(3,2)*b + dc(3,3)*c

      end if

      return
      end


!**********************************************************************
      subroutine inverse (dc,dci)
!**********************************************************************
!**********************************************************************
      double precision ::  dc(3,3),dci(3,3),det,deti
    
      det= dc(1,1)*( dc(2,2)*dc(3,3) - dc(3,2)*dc(2,3) ) &
        - dc(2,1)*( dc(1,2)*dc(3,3) - dc(3,2)*dc(1,3) )  &
        + dc(3,1)*( dc(1,2)*dc(2,3) - dc(2,2)*dc(1,3) )
      deti=1.d0/det

      dci(1,1)= (dc(2,2)*dc(3,3) - dc(3,2)*dc(2,3))*deti
      dci(1,2)=-(dc(1,2)*dc(3,3) - dc(3,2)*dc(1,3))*deti
      dci(1,3)= (dc(1,2)*dc(2,3) - dc(2,2)*dc(1,3))*deti

      dci(2,1)=-(dc(2,1)*dc(3,3) - dc(3,1)*dc(2,3))*deti
      dci(2,2)= (dc(1,1)*dc(3,3) - dc(3,1)*dc(1,3))*deti
      dci(2,3)=-(dc(1,1)*dc(2,3) - dc(2,1)*dc(1,3))*deti

      dci(3,1)= (dc(2,1)*dc(3,2) - dc(3,1)*dc(2,2))*deti
      dci(3,2)=-(dc(1,1)*dc(3,2) - dc(3,1)*dc(1,2))*deti
      dci(3,3)= (dc(1,1)*dc(2,2) - dc(2,1)*dc(1,2))*deti



      return
      end
!**********************************************************************
!**********************************************************************

	     
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) write(*,*)  'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 'k'1k30m,t+W.
	     


