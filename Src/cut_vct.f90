    program cut_vct
       implicit none
       
       real(8), parameter :: a0=2.8553
       integer,dimension(:), allocatable :: n1,n2,n3,nloop,nocto,nloopc,noctoc,ntcube
       double precision,dimension(:),allocatable :: rx,ry,rz,rxb,ryb,rzb
       real(8), dimension(3,3) :: a
       real(8) :: xmin,xmax,ymin,ymax,zmin,zmax,radius
       real(8) :: XMAXC,YMAXC,ZMAXC,XMINC,YMINC,ZMINC
       integer :: tsia,ntemp,nsize,nz16,nz16c,ncube
       integer :: n1MAX,n2MAX,n3MAX,icenter
       
       
       write(*,*) 'Choice what  you want to cut: (1 - 110)(2 - 111)(3 - 100) (4-gao)(5 -vac)(6-octa)(7-photo):'
       read(*,*) tsia
       if (.not.((tsia.eq.6).or.(tsia.eq.7)) )  then 
        write (*,*) 'tsia ......:', tsia, 'not implemented (yet)'
	stop   	
       end if 
       
       if (tsia.eq.6)  then
        !we read the input cell of the perfect bulk ...
        !read the ntemp from formd.str 
	call read_ntemp(ntemp)       
         nsize=ntemp
        allocate (n1(nsize),  n2(nsize),  n3(nsize),        &
                  rx(nsize),  ry(nsize), rz(nsize),         &      
                  rxb(nsize),  ryb(nsize), rzb(nsize),      &
		  nloop(nsize),nocto(nsize),ntcube(nsize),                &
                  nloopc(nsize),noctoc(nsize))           
     
        call read_input(ntemp,nsize,a,rx,ry,rz)
        call write_n123 (ntemp,nsize,a0,n1,n2,n3,rx,ry,rz,  &
                               n1MAX,n2MAX,n3MAX)  
       
	call find_the_center (rx,ry,rz,ntemp,xmin,xmax,ymin,ymax,zmin,zmax,icenter) 
		  
	  write(*,'("The central atoms is ................:",i7)') icenter
	  write(*,'("The coordinates of the central atom..:",3i7)')   n1(icenter),n2(icenter),n3(icenter) 
	  write(*,'("The radius must be smaller than......:",3f12.4)')     &
	              xmax-xmin,ymax-ymin,zmax-zmin
!	  write(*,*) 'Without interaction smaller than.....:'
	  write(*,*) 'Put the radius.......................>'
	  read(*,*) radius
	  
!	  do i=1,ntemp
!	    if ((n1(i)==((n1MAX+1)/2)).and.(n1(i)==((n3MAX+1)/2)).and.(n3(i)==((n3MAX+1)/2)) ) then
!	    icenter=i
!	    end if
!	  end do
  
!	  write(*,*) icenter,n1(icenter),n2(icenter),n3(icenter)
! nloopc, noctoc, nz16c - are for the cluster
! nloop ,  nocto, nz16  - nocto, nloop, ntcube
          call cut_the_bcc (rx,ry,rz,ntemp,nloopc,icenter,radius)
          call find_the_nocto (ntemp,noctoc,nloopc,n1,n2,n3,nz16c)
          call write_no00inp (ntemp,n1,n2,n3,nloopc,noctoc,nz16c) 
          call write_no00xyz (ntemp,rx,ry,rz,n1,n2,n3,nloopc,noctoc,nz16c, &
                                  XMAXC,YMAXC,ZMAXC,XMINC,YMINC,ZMINC )
          nloop(:)=1 
          call find_the_nocto (ntemp,nocto,nloop,n1,n2,n3,nz16)
          call find_the_ncube (ntemp,ntcube,n1,n2,n3,ncube)
          call write_n00_and_beyondxyz(ntemp,rx,ry,rz,n1,n2,n3,nloop,nocto,&
                          nloopc,noctoc,ntcube,nz16,nz16c,ncube , &
                           XMAXC,YMAXC,ZMAXC,XMINC,YMINC,ZMINC,a0 )
          

       end if

       if (tsia.eq.7) then
        call read_ntemp(ntemp)       
         nsize=ntemp
        allocate (n1(nsize),  n2(nsize),  n3(nsize),           &
                  rx(nsize),  ry(nsize), rz(nsize),            &      
                  rxb(nsize),  ryb(nsize), rzb(nsize),         &
                  ntcube(ntemp),                               &
                  nloop(nsize),nocto(nsize))           
     
        call read_input(ntemp,nsize,a,rx,ry,rz)
        call write_n123 (ntemp,nsize,a0,n1,n2,n3,rx,ry,rz,  &
                               n1MAX,n2MAX,n3MAX)  

	  nloop(:)=1
	  call find_the_nocto (ntemp,nocto,nloop,n1,n2,n3,nz16)
	  call find_the_ncube (ntemp,ntcube,n1,n2,n3,ncube)
	  call write_xyz_nocto (ntemp,rx,ry,rz,nloop,nocto,ntcube,nz16,ncube)   
        
       end if

        end  program cut_vct
	
!#####################################################################
!#####################################################################
!#####################################################################
       	  
       subroutine find_the_center (rx,ry,rz,ntemp,xmin,xmax,ymin,ymax,zmin,zmax,icenter)  
       implicit none
        integer, intent(in) ::   ntemp
        real(8), intent(in) ::  rx(ntemp),ry(ntemp),rz(ntemp)
	real(8), intent(out) :: xmin,xmax,ymin,ymax,zmin,zmax
        integer, intent(out) :: icenter
	integer :: i
	real(8) :: xi,yi,zi,ANLIMIT,d(ntemp) 
 
   	  xmax=MAXVAL(rx)
	  xmin=MINVAL(rx)
	  ymax=MAXVAL(ry)
	  ymin=MINVAL(ry)
	  zmax=MAXVAL(rz)
	  zmin=MINVAL(rz)

	  	  
	  ANLIMIT=9999999.0
	  do i=1,ntemp
	   xi=rx(i)-(xmax-xmin)/2.0
	   yi=ry(i)-(ymax-ymin)/2.0
	   zi=rz(i)-(zmax-zmin)/2.0
	   
	    d(i)=xi**2+yi**2+zi**2
	    
	    if (d(i) .lt. ANLIMIT) then
	      icenter=i
	      ANLIMIT=d(i)
	    end if
	  
	  end do
      return	  
      end
!#####################################################################
       	  
	subroutine cut_the_bcc (rx,ry,rz,ntemp,nloop,icenter,radius)
        ! nloop(i) the number of initialy selected atoms from the BCC lattice
	implicit none 
        integer, intent(in) ::   ntemp
        real(8), intent(in) ::  radius,rx(ntemp),ry(ntemp),rz(ntemp)
	integer, intent(in) ::  icenter
	real(8)             :: dtest,radius2,xc,yc,zc
	integer, intent(out):: nloop(ntemp)
	integer :: i, icnt
	
	
	xc=rx(icenter)
	yc=ry(icenter)
	zc=rz(icenter)
	radius2=radius**2
	
	 nloop(:)=0
	 icnt=0
	 do i=1,ntemp
	      dtest=(rx(i)-xc)**2 + (ry(i)-yc)**2 + (rz(i)-zc)**2
	      if (dtest.lt.radius2) then
	       nloop(i)=1
	       icnt=icnt+1
	      end if 
	 end do
	 
	 write(*,'("No of atoms selected initialy from the BCC lattice .....",i7)') icnt
	
	
	
	return
	end
!#####################################################################
       	  
	subroutine find_the_nocto (ntemp,nocto,nloop,n1,n2,n3,nz16)
        ! nocto(i)=1 - the Z16 with  0 rotation
        ! nocto(i)=2 - the Z16 with 90 orientation 
	implicit none 
        integer, intent(in) ::   ntemp
	integer, intent(in):: nloop(ntemp),n1(ntemp),n2(ntemp),n3(ntemp)
	integer, intent(out):: nocto(ntemp),nz16
	integer :: i,icnt1,icnt2
	real(8) :: n1t,n2t,n3t
	
	nocto(:)=0
	icnt1=0
	icnt2=0
	do i=1,ntemp
	 if (nloop(i).eq.1) then
!	   n1t=dble( n1(i)+n2(i)-n3(i))/4.d0
!	   n2t=dble( n1(i)-n2(i)+n3(i))/4.d0
!	   n3t=dble(-n1(i)+n2(i)+n3(i))/4.d0

	   n1t=dble( n1(i)-n2(i)+n3(i)-1)/4.d0
	   n2t=dble( n1(i)+n2(i)-n3(i)-1)/4.d0
	   n3t=dble(-n1(i)+n2(i)+n3(i)-1)/4.d0

	   if ( (dabs(n1t-nint(n1t)).lt.1d-6).and. &
	        (dabs(n2t-nint(n2t)).lt.1d-6).and. &
	        (dabs(n3t-nint(n3t)).lt.1d-6) ) then
	     nocto(i) = 1
	     icnt1=icnt1+1
!debug	     write(*,*) 'type 1', n1(i),n2(i),n3(i)
	   end if
	   if (nocto(i).ne.1) then
	  
	    n1t=n1t-1.d0/4.d0
	    n2t=n2t-1.d0/4.d0
	    n3t=n3t-1.d0/4.d0
	    if ( (dabs(n1t-nint(n1t)).lt.1d-6).and. &
	         (dabs(n2t-nint(n2t)).lt.1d-6).and. &
	         (dabs(n3t-nint(n3t)).lt.1d-6) ) then
	         nocto(i) = 2
	     icnt2=icnt2+1
!debug	       write(*,*) 'type 2', n1(i),n2(i),n3(i)
	    end if
	   end if 
	 
	      
	      
	 end if
	end do
	
      write(*,'("Type 1 and 2  Z16........:",i5)') icnt1,icnt2 	
       nz16=icnt1+icnt2
       
    return
    end	
	
!#####################################################################
       	  
	subroutine find_the_ncube (ntemp,ntcube,n1,n2,n3,ncube) 
        ! ntcube - gives all the sites which are not the center of the cube
	implicit none 
        integer, intent(in) :: ntemp
	integer, intent(in) :: n1(ntemp),n2(ntemp),n3(ntemp)
	integer, intent(out):: ntcube(ntemp),ncube
	integer :: i,icnt
	real(8) :: n1t,n2t,n3t
	
	ntcube(:)=0
	icnt=0
	do i=1,ntemp
	   n1t=dble( n1(i))/2.d0
	   n2t=dble( n2(i))/2.d0
	   n3t=dble( n3(i))/2.d0
	   if ( (dabs(n1t-nint(n1t)).lt.1d-6).and. &
	        (dabs(n2t-nint(n2t)).lt.1d-6).and. &
	        (dabs(n3t-nint(n3t)).lt.1d-6) ) then
	     ntcube(i) = 1
	     icnt=icnt+1
	 
	      
	      
	 end if
	end do
	
      write(*,'("Type Cube.............:",i5)') icnt	
       ncube=icnt       
    return
    end	
	
	
	
!#####################################################################
       	  
	subroutine write_no00inp (ntemp,n1,n2,n3,nloop,nocto,nz16) 
	implicit none
	integer, intent(in) :: ntemp,nz16
	integer, intent(in) :: n1(ntemp),n2(ntemp),n3(ntemp),nloop(ntemp),nocto(ntemp)
	integer   :: i
	
	open (unit=11, file='no00.inp_new', status='unknown')
	write(11,'(i7)') nz16
	
	do i=1,ntemp
	 if ( (nloop(i).eq.1).and.(nocto(i).ne.0) ) then
         write(11,'(3i7)') n1(i),n2(i),n3(i)	
         end if         	
	end do
	
	do i=1,ntemp
	 if ( (nloop(i).eq.1).and.(nocto(i).ne.0) ) then
         write(11,'(3i7)') nocto(i)	
         end if         	
	end do
	
	close(11, status='keep')
	
	return
	end
	
!#####################################################################
       	  
	subroutine write_no00xyz (ntemp,rx,ry,rz,n1,n2,n3,nloop,nocto,nz16, &
                                  XMAX,YMAX,ZMAX,XMIN,YMIN,ZMIN ) 
	implicit none
	integer, intent(in) :: ntemp,nz16
	integer, intent(in) :: nloop(ntemp),nocto(ntemp)
        integer, intent(in) :: n1(ntemp),n2(ntemp),n3(ntemp)
        real(8), intent(in) :: rx(ntemp),ry(ntemp),rz(ntemp)
        real(8), intent(out):: XMAX,YMAX,ZMAX,XMIN,YMIN,ZMIN
	integer   :: i
	
        XMAX=-1d30
        YMAX=-1d30
        ZMAX=-1d30 
        XMIN= 1d30
        YMIN= 1d30
        ZMIN= 1d30 


	open (unit=11, file='center.xyz', status='unknown')
	write(11,'(i7)') nz16
        write(11,'(a)') ' '

	do i=1,ntemp
	 if ( (nloop(i).eq.1).and.(nocto(i).ne.0) ) then
          write (11, '("Ni      ",3f12.5,3i5)') rx(i),ry(i),rz(i), n1(i),n2(i),n3(i)

          if (rx(i).ge.XMAX) XMAX=rx(i)
          if (ry(i).ge.YMAX) YMAX=ry(i)
          if (rz(i).ge.XMAX) ZMAX=rz(i)
          if (rx(i).le.XMIN) XMIN=rx(i)
          if (ry(i).le.YMIN) YMIN=ry(i)
          if (rz(i).le.ZMIN) ZMIN=rz(i)
         end if         	
	end do
        
	
	close(11, status='keep')


	
	return
	end


!#####################################################################
       	  
       subroutine write_xyz_nocto (ntemp,rx,ry,rz,nloop,nocto,ntcube,nz16,ncube) 
       
       implicit none 
       
       integer, intent(in)  :: ntemp,nz16,ncube
       integer,intent(in) ::  nloop(ntemp),nocto(ntemp),ntcube(ntemp)
       double precision, intent(in)     :: rx(ntemp),ry(ntemp),rz(ntemp)
       integer :: i,icnt
       
       open(11,file='jsia.xyz')


       write(11,'(i7)') nz16+ncube+nz16/2
       write(11,'(a)') ' '
       open(12,file='cube.xyz')


       write(12,'(i7)') ntemp
       write(12,'(a)') ' '

       
       icnt=0
	do i=1,ntemp
	
	 if ( (nloop(i).eq.1).and.(nocto(i).ne.0) ) then
	 icnt=icnt+1
           if (nocto(i).eq.1) write (11, '("Zr  ",3f12.5)') rx(i),ry(i),rz(i)       
           if (nocto(i).eq.2) write (11, '("Zr  ",3f12.5)') rx(i),ry(i),rz(i)       
           if (nocto(i).eq.2) write (11, '("Fe  ",3f12.5)') rx(i),ry(i),rz(i)       
         end if
	   if (ntcube(i).eq.1) write(11,'("Cu  ",3f12.5)') rx(i),ry(i),rz(i)               	
	   if (ntcube(i).eq.1) write(12,'("Cu  ",3f12.5)') rx(i),ry(i),rz(i)               	
	   if (ntcube(i).eq.0) write(12,'("Ni  ",3f12.5)') rx(i),ry(i),rz(i)               	
	end do
	
	
	
       close(11)
       


       return
       end
!#####################################################################
	
       subroutine write_n00_and_beyondxyz(ntemp,rx,ry,rz,n1,n2,n3,nloop,nocto,&
                                      nloopc,noctoc,ntcube,nz16,nz16c,ncube , &
                                  XMAXC,YMAXC,ZMAXC,XMINC,YMINC,ZMINC,a0 )
       implicit none 
       
       integer, intent(in)  :: ntemp,nz16,nz16c,ncube
       integer,intent(in) ::  nloop(ntemp),nocto(ntemp),ntcube(ntemp) 
       integer,intent(in) ::  nloopc(ntemp),noctoc(ntemp)
       double precision, intent(in)     :: rx(ntemp),ry(ntemp),rz(ntemp)
       integer , intent(in) :: n1(ntemp),n2(ntemp),n3(ntemp)
       real(8), intent(inout) ::  XMAXC,YMAXC,ZMAXC,XMINC,YMINC,ZMINC
       real(8), intent(in)    :: a0 
       integer :: i,icnt
       real(8)   :: factor
       
       factor=1.1d0
       XMAXC=XMAXC+factor*a0
       YMAXC=YMAXC+factor*a0
       ZMAXC=ZMAXC+factor*a0
       XMINC=XMINC-factor*a0
       YMINC=YMINC-factor*a0
       ZMINC=ZMINC-factor*a0
       icnt=0
        do i=1,ntemp
	 if ( ((rx(i).lt.XMAXC).and.(rx(i).gt.XMINC)).and. &
              ((ry(i).lt.YMAXC).and.(ry(i).gt.YMINC)).and. &
              ((rz(i).lt.ZMAXC).and.(rz(i).gt.ZMINC)) ) then


	  if ( (nloop(i).eq.1).and.(nocto(i).ne.0) ) then
           if (nocto(i).eq.1)  icnt=icnt+1
           if (nocto(i).eq.2)  icnt=icnt+1
           if (nocto(i).eq.2)  icnt=icnt+1 
         end if
          if (ntcube(i).eq.1)  icnt=icnt+1
         end if
        end do 

        do i=1,ntemp
	 if ( (nloopc(i).eq.1).and.(noctoc(i).ne.0) ) then
           icnt=icnt+1
         end if  
        end do
      

       
       open(11,file='center_plus_bulk.xyz')

       
       write(11,'(i7)') icnt
       write(11,'(a)') ' '
       
       icnt=0
        do i=1,ntemp
	 if ( ((rx(i).lt.XMAXC).and.(rx(i).gt.XMINC)).and. &
              ((ry(i).lt.YMAXC).and.(ry(i).gt.YMINC)).and. &
              ((rz(i).lt.ZMAXC).and.(rz(i).gt.ZMINC)) ) then


	  if ( (nloop(i).eq.1).and.(nocto(i).ne.0) ) then
	   icnt=icnt+1
            if (nocto(i).eq.1) write (11, '("Zr  ",3f12.5,3i5)') rx(i),ry(i),rz(i), n1(i),n2(i),n3(i)
       
            if (nocto(i).eq.2) write (11, '("Zr  ",3f12.5,3i5)') rx(i),ry(i),rz(i), n1(i),n2(i),n3(i)
       
            if (nocto(i).eq.2) write (11, '("Fe  ",3f12.5,3i5)') rx(i),ry(i),rz(i), n1(i),n2(i),n3(i)
       
           end if
	    if (ntcube(i).eq.1) write(11,'("Cu  ",3f12.5,3i5)') rx(i),ry(i),rz(i), n1(i),n2(i),n3(i)
               	
         end if
        end do

	do i=1,ntemp
	 if ( (nloopc(i).eq.1).and.(noctoc(i).ne.0) ) then
          write (11, '("Ni      ",3f12.5,3i5)') rx(i),ry(i),rz(i), n1(i),n2(i),n3(i)
         end if         	
	end do
       
       close(11)
       


       return
       end

