

       implicit double precision (a-h,o-z)
       parameter (a0=2.8553,as110=2.148d0,as111=2.12d0,as100=2.1d0)
       integer,dimension(:), allocatable :: n1(:),n2(:),n3(:)
       double precision,dimension(:),allocatable :: rx(:),ry(:),rz(:),&
                                                 rxc(:),ryc(:),rzc(:),&
					rxloop(:),ryloop(:),rzloop(:),&
			             rxloop2(:),ryloop2(:),rzloop2(:)						 
       integer,dimension(:),allocatable :: ns1(:),ns2(:),ns3(:),nsite(:),&
                                           nicnt  
       double precision     :: rsx(2),rsy(2),rsz(2)
      double precision     :: nv(3),                                     &
                              nv1(3),nv2(3),nv3(3),nv4(3),nv5(3),nv6(3)

       double precision,dimension(:,:) :: a(3,3)
       integer :: tsia,nsia,hplane,bshape,biplane
       
       
!       write(*,*) 'Choice the sia (1 - 110)(2 - 111)(3 - 100):'
!       read(*,*) tsia
!       write(*,*) 'excelent choice master'
       
       call read_ntemp(nsolid)
       
       print*,  'Total number of the perfect solid....:',nsolid

       
       allocate (rx(nsolid)    ,  ry(nsolid)    , rz(nsolid)    ,   &
                 rxc(nsolid)   ,  ryc(nsolid)   , rzc(nsolid)   ,   &
                 rxloop(nsolid),  ryloop(nsolid), rzloop(nsolid),   &
                 nicnt(nsolid) )

       call read_input(nsolid,a,rx,ry,rz)
       
       write(*,*) 'The size of the box..................:'
       write(*,'(3f12.5)')  a   

       
       write(*,*) 'If the loop is in a biplane or plane:'
       write(*,*) '1/plane, 2/bi-plane..................> '
       read(*,*)  biplane
       write(*,*) 'The habit plane......................:'
       write(*,*) '1 - 100/ 2 - 111 / 3  - 110..........>'
       read (*,*) hplane 
       
       if (hplane==1) then
!        normal of the plane 

	   nv(1)= 0.0
	   nv(2)= 0.0
	   nv(3)= 1.0

        call cut_plane(nsolid,nplane,nv,nicnt,rx,ry,rz,rxc,ryc,rzc)

       write(*,*) 'The shape of the loop................:'
       write(*,*) '(1/rectangle, 2/circle, 3/hexa)......>'
         read(*,*) bshape

	 if (bshape==1) then
	  rxmax=MAXVAL(rx)/2.0
	  rymax=MAXVAL(ry)/2.0

	  write(*,*) 'The X side must be smaller than......:',rxmax
	  write(*,*) 'X side...............................>'
          read(*,*) side1
	  write(*,*) 'The Y side must be smaller than......:',rymax
	  write(*,*) 'Y side...............................>'
          read(*,*) side2
	   
	   
	   
	   
	   nv1(1)= 1.0
	   nv1(2)= 0.0
	   nv1(3)= 0.0

	 
	 call count_cut_rec(nsolid,side1,side2,nplane,iloop,         &
	              nv,nv1,&
	              nicnt,&
                      rxc,ryc,rzc,                                &
		      rxloop,ryloop,rzloop)       
         end if 
	 
	 if (bshape==2) then
	  rmax=MAXVAL(rx)/2.0
	  write(*,*) 'The radius must be smaller than......:',rmax
	  write(*,*) 'Without interaction smaller than.....:',rmax/2.0
	  write(*,*) 'Put the radius.......................>'
          read(*,*) rcirc
	 call count_cut_circ(nsolid,rcirc,nplane,iloop,nicnt,&
                      rxc,ryc,rzc,                           &
		      rxloop,ryloop,rzloop)       
         end if
	 
	 if (bshape==3) then

	  xmax=MAXVAL(rx)/2.0
	  xmin=MINVAL(rx)/2.0
	  ymax=MAXVAL(ry)/2.0
	  ymin=MINVAL(ry)/2.0
	  zmax=MAXVAL(rz)/2.0
	  zmin=MINVAL(rz)/2.0
	  write(*,*) 'The radius must be smaller than......:',     &
	              xmax-xmin,ymax-ymin,zmax-zmin
!	  write(*,*) 'Without interaction smaller than.....:'
	  write(*,*) 'Put the radius.......................>'
          read(*,*) side
	         write(*,*) nv	
	nv1(1)=  1.0
	nv1(2)=  0.0
	nv1(3)=  0.0

	nv2(1)=  0.0
	nv2(2)=  1.0
	nv2(3)=  0.0

	nv3(1)=  1.0/sqrt(2.0)
	nv3(2)=  1.0/sqrt(2.0)
	nv3(3)=  0.0

	nv4(1)= -1.0/sqrt(2.0) 
	nv4(2)=  1.0/sqrt(2.0) 
	nv4(3)=  0.0

	nv5(1)=  0.0
	nv5(2)=  0.0
	nv5(3)=  1.0

	nv6(1)=  0.0
	nv6(2)=  0.0
	nv6(3)= -1.0

	 
	 
	 call count_cut_hex(nsolid,side,nplane,iloop,      &
	 nv,nv1,nv2,nv3,nv4,nv5,nv6,                           &
	              nicnt,&
                      rxc,ryc,rzc,                             &
		      rxloop,ryloop,rzloop)
        end if



       end if  
       

       if (hplane==2) then
!        normal of the plane 

        nv(1)=1.d0/dsqrt(3.d0) 
        nv(2)=1.d0/dsqrt(3.d0) 
        nv(3)=1.d0/dsqrt(3.d0) 


        call cut_plane(nsolid,nplane,nv,nicnt,rx,ry,rz,rxc,ryc,rzc)

         write(*,*) 'The shape of the loop................:'
         write(*,*) '(1/rectangle, 2/circle, 3/hexa)......>'
         read(*,*) bshape
	 
	 
	if (bshape==1) then
	  rxmax=MAXVAL(rx)/2.0
	  rymax=MAXVAL(ry)/2.0

	  write(*,*) 'The X side must be smaller than......:',rxmax
	  write(*,*) 'X side...............................>'
          read(*,*) side1
	  write(*,*) 'The Y side must be smaller than......:',rymax
	  write(*,*) 'Y side...............................>'
          read(*,*) side2
  	   
	   nv1(1)= 1.0/sqrt(2.0)
	   nv1(2)=-1.0/sqrt(2.0)
	   nv1(3)= 0.0
	  
	 call count_cut_rec(nsolid,side1,side2,nplane,iloop,nv,nv1,  &
	              nicnt,&
                      rxc,ryc,rzc,                                &
		      rxloop,ryloop,rzloop)       
         end if 

	 if (bshape==2) then
	  xmax=MAXVAL(rx)/2.0
	  xmin=MINVAL(rx)/2.0
	  ymax=MAXVAL(ry)/2.0
	  ymin=MINVAL(ry)/2.0
	  zmax=MAXVAL(rz)/2.0
	  zmin=MINVAL(rz)/2.0
	  write(*,*) 'The radius must be smaller than......:',     &
	              xmax-xmin,ymax-ymin,zmax,zmin
!	  write(*,*) 'Without interaction smaller than.....:'
	  write(*,*) 'Put the radius.......................>'
          read(*,*) rcirc
	 call count_cut_circ(nsolid,rcirc,nplane,iloop,nicnt,&
                      rxc,ryc,rzc,                           &
		      rxloop,ryloop,rzloop)       
         end if
	 
	 if (bshape==3) then

	  xmax=MAXVAL(rx)/2.0
	  xmin=MAXVAL(rx)/2.0
	  ymax=MAXVAL(ry)/2.0
	  ymin=MAXVAL(ry)/2.0
	  zmax=MAXVAL(rz)/2.0
	  zmin=MAXVAL(rz)/2.0
	  write(*,*) 'The radius must be smaller than......:',     &
	              xmax-xmin,ymax-ymin,zmax,zmin
!	  write(*,*) 'Without interaction smaller than.....:'
	  write(*,*) 'Put the radius.......................>'
          read(*,*) side
	 
          	
	 nv1(1)=  1.0
	 nv1(2)=  0.0
	 nv1(3)=  0.0

	 nv2(1)= -1.0
	 nv2(2)=  0.0
	 nv2(3)=  0.0

	 nv3(1)=  0.0
	 nv3(2)=  1.0
	 nv3(3)=  0.0

	 nv4(1)=  0.0
	 nv4(2)= -1.0
	 nv4(3)=  0.0

	 nv5(1)=  0.0
	 nv5(2)=  0.0
	 nv5(3)=  1.0

	 nv6(1)=  0.0
	 nv6(2)=  0.0
	 nv6(3)= -1.0
	 
	 
	 call count_cut_hex(nsolid,side,nplane,iloop,      &
	 nv,nv1,nv2,nv3,nv4,nv5,nv6,                       &
	              nicnt,&
                      rxc,ryc,rzc,                         &
		      rxloop,ryloop,rzloop)

         end if


       end if 
       
       if (hplane==3) then
	   nv(1)= 1.0/sqrt(2.0)
	   nv(2)= 1.0/sqrt(2.0)
	   nv(3)= 0.0
	
	call cut_plane(nsolid,nplane,nv,nicnt,rx,ry,rz,rxc,ryc,rzc)       
	   
         write(*,*) 'The shape of the loop................:'
         write(*,*) '(1/rectangle, 2/circle, 3/hexa)......>'
         read(*,*) bshape

        if (bshape==1) then
	 rxmax=MAXVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	 rxmin=MINVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	 rymax=MAXVAL(rz(:))/2.0
	 rymin=MINVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	 rxmin=MINVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	  write(*,*) 'The Z side must be smaller than......:',rymax-rymin
	  write(*,*) 'Z side...............................>'
          read(*,*) side2
	  write(*,*) 'The XY side must be smaller than.....:',rxmax-rxmin
	  write(*,*) 'X side...............................>'
          read(*,*) side1
	
	   nv1(1)= 1.0
	   nv1(2)= 0.0
	   nv1(3)= 0.0
	  
	 call count_cut_rec(nsolid,side1,side2,nplane,iloop,nv,nv1,  &
	              nicnt,&
                      rxc,ryc,rzc,                                &
		      rxloop,ryloop,rzloop)       
	end if
	 
	 if (bshape==2) then
	 rxmax=MAXVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	 rxmin=MINVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	 rymax=MAXVAL(rz(:))/2.0
	 rxmin=MINVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	 rymin=MINVAL(sqrt(rx(:)**2+ry(:)**2))/2.0

	  write(*,*) 'The radius must be smaller than......:',     &
	              rxmax-rxmin,rymax-rymin
!	  write(*,*) 'Without interaction smaller than.....:'
	  write(*,*) 'Put the radius.......................>'
          read(*,*) rcirc
	  call count_cut_circ(nsolid,rcirc,nplane,iloop,nicnt,&
                      rxc,ryc,rzc,                           &
		      rxloop,ryloop,rzloop)       
         end if

	
	 if (bshape==3) then

	 rxmax=MAXVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	 rxmin=MINVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	 rymax=MAXVAL(rz(:))/2.0
	 rymin=0
	 rxmin=MINVAL(sqrt(rx(:)**2+ry(:)**2))/2.0
	  write(*,*) 'The radius must be smaller than......:',     &
	              rxmax-rxmin,rymax-rymin
!	  write(*,*) 'Without interaction smaller than.....:'
	  write(*,*) 'Put the radius.......................>'
          read(*,*) side
	 
          	
!	 nv1(1)= -1.0/sqrt(2.0)
!	 nv1(2)=  1.0/sqrt(2.0)
!	 nv1(3)=  0.0

!	 nv2(1)=  1.0/sqrt(2.0)
!	 nv2(2)= -1.0/sqrt(2.0)
!	 nv2(3)=  0.0
	 
!	 ca2=(1.d0+1.d0/dsqrt(3.d0))/2.d0
!	 ca=dsqrt(ca2)
!	 sa=dsqrt(1.d0-ca2)
!	 ana=1.d0/dsqrt(2.d0)

!	 nv3(1)=  sa*ana
!	 nv3(2)= -sa*ana
!	 nv3(3)=  ca
	 
!	 print *, nv3(1)**2+nv3(2)**2+nv3(3)**2

!	 nv4(1)= -sa*ana
!	 nv4(2)=  sa*ana
!	 nv4(3)= -ca

!	 nv5(1)= -sa*ana
!	 nv5(2)=  sa*ana
!	 nv5(3)=  ca

!	 nv6(1)=  sa*ana
!	 nv6(2)= -sa*ana
!	 nv6(3)= -ca


	nv1(1)= -1.0/sqrt(2.0)
	nv1(2)=  1.0/sqrt(2.0)
	nv1(3)=  0.0

	nv2(1)=  1.0/sqrt(2.0)
	nv2(2)= -1.0/sqrt(2.0)
	nv2(3)=  0.0
        
        ab=1.d0/sqrt(6.d0)
	cc=sqrt(2.d0/3.d0)
	
	nv3(1)= -ab
	nv3(2)=  ab
	nv3(3)=  cc
        

	nv4(1)=  ab
	nv4(2)= -ab
	nv4(3)= -cc

	nv5(1)=  ab
	nv5(2)= -ab
	nv5(3)=  cc

	nv6(1)= -ab
	nv6(2)=  ab
	nv6(3)= -cc
	 
	 
	 call count_cut_hex2(nsolid,side,nplane,iloop,      &
	 nv,nv1,nv2,nv3,nv4,nv5,nv6,                       &
	              nicnt,&
                      rxc,ryc,rzc,                         &
		      rxloop,ryloop,rzloop)

         end if

	end if
       




       

!.....................................................................
       if (biplane==2) then
            nsia=2*iloop
            nall=nsolid+nsia
 	    allocate (n1(nsolid),  n2(nsolid),  n3(nsolid), &
 		  ns1(nsia),ns2(nsia),ns3(nsia),nsite(nsia), &
 		  rxloop2(nsia),ryloop2(nsia),rzloop2(nsia) )  

 	    if (hplane==1) then

 	     do i=1,nsia/2
 	     rxloop2(i)= rxloop(i)
 	     ryloop2(i)= ryloop(i)
 	     rzloop2(i)= rzloop(i)
 	     end do
 	     do i=nsia/2+1,nsia
 	     rxloop2(i)= rxloop(i-nsia/2)+a0/2
 	     ryloop2(i)= ryloop(i-nsia/2)+a0/2
 	     rzloop2(i)= rzloop(i-nsia/2)+a0/2       
 	     end do
	    else if (hplane==2) then
 	     do i=1,nsia/2
 	     rxloop2(i)= rxloop(i)
 	     ryloop2(i)= ryloop(i)
 	     rzloop2(i)= rzloop(i)
 	     end do
 	     do i=nsia/2+1,nsia
 	     rxloop2(i)= rxloop(i-nsia/2)+a0/2
 	     ryloop2(i)= ryloop(i-nsia/2)+a0/2
 	     rzloop2(i)= rzloop(i-nsia/2)-a0/2       
 	     end do
	    
            else if (hplane==3) then
 	     do i=1,nsia/2
 	     rxloop2(i)= rxloop(i)
 	     ryloop2(i)= ryloop(i)
 	     rzloop2(i)= rzloop(i)
 	     end do
 	     do i=nsia/2+1,nsia
 	     rxloop2(i)= rxloop(i-nsia/2)+a0/2
 	     ryloop2(i)= ryloop(i-nsia/2)
 	     rzloop2(i)= rzloop(i-nsia/2)+a0/2       
 	     end do
	     
 	    
 	    end if
       write(*,*)    'The size of the loop (no of atoms)...:' , nsia

         else if (biplane==1) then
            nsia=iloop
            nall=nsolid+nsia
 	    allocate (n1(nsolid),  n2(nsolid),  n3(nsolid), &
 		  ns1(nsia),ns2(nsia),ns3(nsia),nsite(nsia), &
 		  rxloop2(nsia),ryloop2(nsia),rzloop2(nsia) )  

 	     do i=1,nsia
 	     rxloop2(i)= rxloop(i)
 	     ryloop2(i)= ryloop(i)
 	     rzloop2(i)= rzloop(i)
 	     end do
 	    
       write(*,*)    'The size of the loop (no of atoms)...:' , nsia
 	    end if

!.....................................................................
!.....................................................................
!.....................................................................

       call write_n123 (nsolid,a0,a,n1,n2,n3,rx,ry,rz,  &
                             n1MAX,n2MAX,n3MAX)  
      
       call config_sia (nsolid,nsia,n1MAX,n2MAX,n3MAX,a0,&
                        ns1,ns2,ns3, &
                        rxloop2,ryloop2,rzloop2)

       
       write(*,*) nsolid, nsia
!       write(*,*) nsite
       call check_site (nsolid,nsia,n1,n2,n3,ns1,ns2,ns3,nsite)       
       

       volc= (dble(nall)/dble(nsolid))**0.3333333333d0
       print *,      'The volume constant factor is........:',  volc         

       write(*,*)    'Choice the sia.......................:'
       write(*,*)    '(1/<110>/2/<111>/3<100>..............>'
              read(*,*) tsia
       write(*,*)    '=======excelent choice master========='

       call type_sia (tsia,hplane,as110,as111,as100,rsx,rsy,rsz)
       
       call write_out (nall,nsolid,nsia,volc,nsite,rx,  &
                       ry,rz,rsx,rsy,rsz,a)


!       call read_sia (tsia,as110,as111,as100,rsx,rsy,rsz)
!       call write_n123 (nsolid,ntemp,a0,a,n1,n2,n3,rx,ry,rz,  &
!                             n1MAX,n2MAX,n3MAX)  
!       call read_config_sia (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
!			     
!       call check_site (nsolid,ntemp,nsia,n1,n2,n3,ns1,ns2,ns3,nsite)
!       volc= (dble(nsolid)/dble(ntemp))**0.3333333333d0
!       print *, volc		         
!       call write_out (nsolid,ntemp,nsia,volc,n1,n2,n3,nsite,rx,ry, &
!                       rz,rsx,rsy,rsz,a)
		
		
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
       	  
       subroutine read_input(nsolid,a,rx,ry,rz)      
       integer              :: nsolid
       double precision     :: rx(nsolid),ry(nsolid),rz(nsolid)
       double precision,dimension(:,:) :: a(3,3)
       character*3 :: sym
       
        read(10,'(3f18.13)')a(1,1), a(1,2), a(1,3)
        read(10,'(3f18.13)')a(2,1), a(2,2), a(2,3)
        read(10,'(3f18.13)')a(3,1), a(3,2), a(3,3)
	
	do i=1,nsolid
		read(10,'(a3,3(f18.13))') sym, rx(i),ry(i),rz(i)
        end do 		
     
		
       
       return
       end 
       	  
!#####################################################################
       	  
       subroutine cut_plane(nsolid,nplane,nv,nicnt,rx,ry,rz,rxc,ryc,rzc)      
       implicit none
       integer              :: nsolid, nplane,i,icnt,im
       integer              :: nicnt(nsolid)
       double precision     ::  rx(nsolid),ry(nsolid),rz(nsolid), &
                             rxc(nsolid),ryc(nsolid),rzc(nsolid), &
			      zplane 
       double precision     :: xmax,xmin,ymax,ymin,zmax,zmin,     &
                               xi,yi,zi,                 &
			       dplane,eqplane,dmin,d(nsolid),nv(3)   			      
!....... find the center 
        xmax=MAXVAL(rx)
	xmin=MINVAL(rx)	
        ymax=MAXVAL(ry)
	ymin=MINVAL(ry)
        zmax=MAXVAL(rz)
	zmin=MINVAL(rz)
	


	do i=1,nsolid
	  xi=rx(i)-(xmax-xmin)/2.0
	  yi=ry(i)-(ymax-ymin)/2.0
	  zi=rz(i)-(zmax-zmin)/2.0
	  d(i)=xi**2+yi**2 +zi**2
	end do	
	
	dmin=MINVAL(d)
	

	do i=1,nsolid
         if (dabs(d(i)-dmin).lt.0.01) then
	    im=i
	 end if        
	end do	
	    write(*,'("The center of the CUBE is............:",i7,3f12.5)') &
	                im,rx(im),ry(im),rz(im)
!....the center was found in the atom number im
       

        dplane=-nv(1)*rx(im)-nv(2)*ry(im)-nv(3)*rz(im)
	 
	write(*,*)'Cutted in the plane ................:',       &
	           -dplane/sqrt(nv(1)**2+nv(2)**2+nv(3)**2)  


	icnt=0
	
	do i=1,nsolid
         eqplane=nv(1)*rx(i)+nv(2)*ry(i)+nv(3)*rz(i)+dplane	 
	 if (dabs(eqplane).gt.0.01) then
!	   ni(i)=0
          else
	   icnt=icnt+1
	   
	  
!  	   write(*,'(2i5,8f12.5)') i,icnt,eqplane,nx,ny,nz,rx(i),ry(i),rz(i),dplane
	   rxc(icnt)=rx(i)
	   ryc(icnt)=ry(i)
   	   rzc(icnt)=rz(i)
           nicnt(icnt)=i
	 if (i.eq.im) then
	   write(12,'("Cu   ",3f12.5)') rxc(icnt),ryc(icnt),rzc(icnt)
	   else
	   write(12,'("Fe   ",3f12.5)') rxc(icnt),ryc(icnt),rzc(icnt)
	end if   
	 end if
	end do	

       nplane=icnt
       write(*,*) 'In this plane the number of atoms....:', nplane
       return
       end 
!#####################################################################

!#####################################################################
       subroutine count_cut_rec(nsolid,side1,side2,nplane,iloop,     &
                      nv,nv1,					     &
                      nicnt,					     &
                      rxc,ryc,rzc,                                   &
		      rxloop,ryloop,rzloop)       
       implicit none
       integer              :: nsolid, nplane
       integer              :: i,im,iloop,nicnt(nsolid)
       double precision     :: nv(3),nv1(3),nv2(3)
       
       double precision     :: rx(nsolid),ry(nsolid),rz(nsolid),    &
                               rxc(nsolid),ryc(nsolid),rzc(nsolid), &
			rxloop(nsolid),ryloop(nsolid),rzloop(nsolid),&
			       xi,yi,zi,side1,side2
       double precision   :: xmin,xmax,ymin,ymax,zmax,zmin,dmin,d(nplane)			        
			       
       
!....... find the center 
        xmax=MAXVAL(rxc)
	xmin=MINVAL(rxc)	
        ymax=MAXVAL(ryc)
	ymin=MINVAL(ryc)
        zmax=MAXVAL(ryc)
	zmin=MINVAL(ryc)
	


	do i=1,nplane
	xi=rxc(i)-(xmax-xmin)/2.0
	yi=ryc(i)-(ymax-ymin)/2.0
	zi=rzc(i)-(zmax-zmin)/2.0
	
	
	d(i)=xi**2+yi**2+zi**2
	end do	
	
	dmin=MINVAL(d)
	

	do i=1,nplane
         if (dabs(d(i)-dmin).lt.0.01) then
	    im=i
	 end if        
	end do	
	    write(*,*) 'The center of the loop is............:', im
!....the center was found in the atom number im
	

        call VectProd(nv,nv1,nv2)

	iloop=0
	do i=1,nplane
	 xi=(rxc(i)-rxc(im))*nv1(1)+(ryc(i)-ryc(im))*nv1(2)+(rzc(i)-rzc(im))*nv1(3)
	 yi=(rxc(i)-rxc(im))*nv2(1)+(ryc(i)-ryc(im))*nv2(2)+(rzc(i)-rzc(im))*nv2(3)
	 if ((dabs(xi).lt.side1/2.0).and.(dabs(yi).lt.side2/2.0)) then
          iloop=iloop+1 
!	   ni(i)=0
           rxloop(iloop)=rxc(i)
           ryloop(iloop)=ryc(i)
           rzloop(iloop)=rzc(i)
	   
!!!	   write(12,*) 'Be  ',rxc(i),ryc(i),rzc(i)
	   write(12,*) rxc(i),ryc(i)

	 end if
	
	end do	
       
       return
       end 

!#####################################################################
       subroutine count_cut_hex(nsolid,side,nplane,iloop,      &
	 nv,nv1,nv2,nv3,nv4,nv5,nv6,                           &
	              nicnt,&
                      rxc,ryc,rzc,                             &
		      rxloop,ryloop,rzloop)
       implicit none
       integer              :: nsolid, nplane
       integer              :: i,im,iloop,nicnt(nsolid)
       double precision     :: nv(3),                                 &
                              nv1(3),nv2(3),nv3(3),nv4(3),nv5(3),nv6(3)
       
       double precision     :: rx(nsolid),ry(nsolid),rz(nsolid),    &
                               rxc(nsolid),ryc(nsolid),rzc(nsolid), &
			rxloop(nsolid),ryloop(nsolid),rzloop(nsolid),&
			       x1i,x2i,x3i,x4i,x5i,x6i,side,xi,yi,zi
       double precision   :: xmin,xmax,ymin,ymax,zmax,zmin,dmin,d(nplane)			        
			       
       
!....... find the center 
        xmax=MAXVAL(rxc)
	xmin=MINVAL(rxc)	
        ymax=MAXVAL(ryc)
	ymin=MINVAL(ryc)
        zmax=MAXVAL(ryc)
	zmin=MINVAL(ryc)
	


	do i=1,nplane
	xi=rxc(i)-(xmax-xmin)/2.0
	yi=ryc(i)-(ymax-ymin)/2.0
	zi=rzc(i)-(zmax-zmin)/2.0
	
	
	d(i)=xi**2+yi**2+zi**2
	end do	
	
	dmin=MINVAL(d)
	

	do i=1,nplane
         if (dabs(d(i)-dmin).lt.0.01) then
	    im=i
	 end if        
	end do	
	    write(*,*) 'The center of the loop is............:', im
!....the center was found in the atom number im
	


	iloop=0
	do i=1,nplane
	 x1i=(rxc(i)-rxc(im))*nv1(1)+(ryc(i)-ryc(im))*nv1(2)+(rzc(i)-rzc(im))*nv1(3)
	 x2i=(rxc(i)-rxc(im))*nv2(1)+(ryc(i)-ryc(im))*nv2(2)+(rzc(i)-rzc(im))*nv2(3)
	 x3i=(rxc(i)-rxc(im))*nv3(1)+(ryc(i)-ryc(im))*nv3(2)+(rzc(i)-rzc(im))*nv3(3)
	 x4i=(rxc(i)-rxc(im))*nv4(1)+(ryc(i)-ryc(im))*nv4(2)+(rzc(i)-rzc(im))*nv4(3)
	 x5i=(rxc(i)-rxc(im))*nv5(1)+(ryc(i)-ryc(im))*nv5(2)+(rzc(i)-rzc(im))*nv5(3)
	 x6i=(rxc(i)-rxc(im))*nv6(1)+(ryc(i)-ryc(im))*nv6(2)+(rzc(i)-rzc(im))*nv6(3)
	 if ((dabs(x1i).lt.side/2.0).and.&
	     (dabs(x2i).lt.side/2.0).and.&
	     (dabs(x3i).lt.side/2.0).and.&
	     (dabs(x4i).lt.side/2.0).and.&
	     (dabs(x5i).lt.side/2.0).and.&
	     (dabs(x6i).lt.side/2.0)) then
          iloop=iloop+1 
!	   ni(i)=0
           rxloop(iloop)=rxc(i)
           ryloop(iloop)=ryc(i)
           rzloop(iloop)=rzc(i)
	   
!!!	   write(12,*) 'Be  ',rxc(i),ryc(i),rzc(i)
	   write(12,*) rxc(i),ryc(i)

	 end if
	
	end do	
       
       return
       end 
!#####################################################################
       	  
       subroutine count_cut_hex2(nsolid,side,nplane,iloop,      &
	 nv,nv1,nv2,nv3,nv4,nv5,nv6,                           &
	              nicnt,&
                      rxc,ryc,rzc,                             &
		      rxloop,ryloop,rzloop)
       implicit none
       integer              :: nsolid, nplane
       integer              :: i,im,iloop,nicnt(nsolid)
       double precision     :: nv(3),                                 &
                              nv1(3),nv2(3),nv3(3),nv4(3),nv5(3),nv6(3)
       
       double precision     :: rx(nsolid),ry(nsolid),rz(nsolid),    &
                               rxc(nsolid),ryc(nsolid),rzc(nsolid), &
			rxloop(nsolid),ryloop(nsolid),rzloop(nsolid),&
			       x1i,x2i,x3i,x4i,x5i,x6i,side,xi,yi,zi
       double precision   :: xmin,xmax,ymin,ymax,zmax,zmin,dmin,d(nplane)			        
			       
       
!....... find the center 
        xmax=MAXVAL(rxc)
	xmin=MINVAL(rxc)	
        ymax=MAXVAL(ryc)
	ymin=MINVAL(ryc)
        zmax=MAXVAL(ryc)
	zmin=MINVAL(ryc)
	


	do i=1,nplane
	xi=rxc(i)-(xmax-xmin)/2.0
	yi=ryc(i)-(ymax-ymin)/2.0
	zi=rzc(i)-(zmax-zmin)/2.0
	
	
	d(i)=xi**2+yi**2+zi**2
	end do	
	
	dmin=MINVAL(d)
	

	do i=1,nplane
         if (dabs(d(i)-dmin).lt.0.01) then
	    im=i
	 end if        
	end do	
	    write(*,*) 'The center of the loop is............:', im
!....the center was found in the atom number im
	


	iloop=0
	do i=1,nplane
	 x1i=(rxc(i)-rxc(im))*nv1(1)+(ryc(i)-ryc(im))*nv1(2)+(rzc(i)-rzc(im))*nv1(3)
	 x2i=(rxc(i)-rxc(im))*nv2(1)+(ryc(i)-ryc(im))*nv2(2)+(rzc(i)-rzc(im))*nv2(3)
	 x3i=(rxc(i)-rxc(im))*nv3(1)+(ryc(i)-ryc(im))*nv3(2)+(rzc(i)-rzc(im))*nv3(3)
	 x4i=(rxc(i)-rxc(im))*nv4(1)+(ryc(i)-ryc(im))*nv4(2)+(rzc(i)-rzc(im))*nv4(3)
	 x5i=(rxc(i)-rxc(im))*nv5(1)+(ryc(i)-ryc(im))*nv5(2)+(rzc(i)-rzc(im))*nv5(3)
	 x6i=(rxc(i)-rxc(im))*nv6(1)+(ryc(i)-ryc(im))*nv6(2)+(rzc(i)-rzc(im))*nv6(3)
	 if ((dabs(x1i).lt.side/2.0).and.&
	     (dabs(x2i).lt.side/2.0).and.&
	     (dabs(x3i).lt.side/2.0).and.&
	     (dabs(x4i).lt.side/2.0).and.&
	     (dabs(x5i).lt.side/2.0).and.&
	     (dabs(x6i).lt.side/2.0)) then
!	     if ((dabs(x1i).lt.side/2.0).and.&
!	     (dabs(x2i).lt.side/2.0)) then
          iloop=iloop+1 
!	   ni(i)=0
           rxloop(iloop)=rxc(i)
           ryloop(iloop)=ryc(i)
           rzloop(iloop)=rzc(i)
	   
!!!!	   write(12,*) 'Be  ',rxc(i),ryc(i),rzc(i)
	   write(12,*) rxc(i),ryc(i)
!           end if
	 end if
	
	end do	
       
       return
       end 
!#####################################################################
       subroutine count_cut_circ(nsolid,rcirc,nplane,iloop,nicnt,&
                      rxc,ryc,rzc,                               &
		      rxloop,ryloop,rzloop)         
       integer              :: nsolid, nplane
       integer              :: i,im,iloop,nicnt(nsolid)

       double precision     :: rx(nsolid),ry(nsolid),rz(nsolid),    &
                      rxloop(nsolid),ryloop(nsolid),rzloop(nsolid), &
		            rxc(nsolid),ryc(nsolid),rzc(nsolid),    &
			       d(nplane), &
			       xi,yi,rcirc,ri,xmin,xmax,ymin,ymax,  &
			       zmax,zmin,zi,dmin
       
!....... find the center 
        xmax=MAXVAL(rxc)
	xmin=MINVAL(rxc)	
        ymax=MAXVAL(ryc)
	ymin=MINVAL(ryc)
        zmax=MAXVAL(rzc)
	zmin=MINVAL(rzc)
	


	do i=1,nplane
	xi=rxc(i)-(xmax-xmin)/2.0
	yi=ryc(i)-(ymax-ymin)/2.0
	zi=rzc(i)-(zmax-zmin)/2.0
	d(i)=xi**2+yi**2+zi**2
	end do	
	
	dmin=MINVAL(d)
	

	do i=1,nplane
         if (dabs(d(i)-dmin).lt.0.01) then
	    im=i
	 end if        
	end do	
	    write(*,*) 'The center of the loop is............:', im
!....the center was found in the atom number im
	
	iloop=0
	do i=1,nplane
	 xi=rxc(i)-rxc(im)
	 yi=ryc(i)-ryc(im)
	 zi=rzc(i)-rzc(im)
	 ri=sqrt(xi**2+yi**2+zi**2)
	 if (ri.lt.rcirc) then
          iloop=iloop+1 
!	   ni(i)=0
           rxloop(iloop)=rxc(i)
           ryloop(iloop)=ryc(i)
           rzloop(iloop)=rzc(i)
!!!	   write(12,*) 'Be  ',rxc(i),ryc(i),rzc(i)
	   write(12,*) rxc(i),ryc(i)

	 end if
	
	end do	
       
       return
       end 
!#####################################################################
       	  
       subroutine type_sia(tsia,hplane,as110,as111,as100,rsx,rsy,rsz)      
       double precision     :: rsx(2),rsy(2),rsz(2),as110,as111,as100
       integer :: tsia,hplane,bshape
       
       cos45 = sqrt(2.d0)/2.d0
       
       
       if (tsia==1) then 
          as2 = as110/2.d0
          rsx(1)= - cos45 * as2
	  rsy(1)= - cos45 * as2
	  rsz(1)= 0.D0
	  
	  rsx(2)= cos45 * as2
	  rsy(2)= cos45 * as2
	  rsz(2)= 0.d0
	   if (hplane==1) then
	    rsx(1)= 0.D0
	    rsy(1)= - cos45 * as2
            rsz(1)= - cos45 * as2
	  
	    rsx(2)= 0.d0
	    rsy(2)= cos45 * as2
	    rsz(2)= cos45 * as2
	   
	   end if
	 else if (tsia==2) then
          as2 = as111/2.d0
          rsx(1)= - cos45 * cos45 *  as2
	  rsy(1)= - cos45 * cos45 *  as2
	  rsz(1)= - cos45 * cos45 *  as2
	  
	  rsx(2)= cos45 * cos45 *  as2
	  rsy(2)= cos45 * cos45 *  as2
	  rsz(2)= cos45 * cos45 *  as2


	 else if (tsia==3) then
          rsx(1)= 0.d0
	  rsy(1)= 0.d0
	  rsz(1)= -as100/2.d0
	  
	  rsx(2)= 0.d0
	  rsy(2)= 0.d0
	  rsz(2)= as100/2.d0
	  if (hplane==3) then
           rsx(1)=-as100/2.d0
	   rsy(1)= 0.d0
	   rsz(1)= 0.d0
	  
	   rsx(2)= as100/2.d0
	   rsy(2)= 0.d0
	   rsz(2)= 0.d0
	  end if 
	 else
	   print*, 'Unknown type of SIA: Check read_sia' 
	   stop
	 end if    
		
        
       return
       end 

!#####################################################################
       	  
       
       
       subroutine write_n123(nsolid,a0,a,n1,n2,n3,rx,ry,rz,  &
                             n1MAX,n2MAX,n3MAX)      
       implicit none
       integer  :: nsolid
       integer  :: n1(nsolid),n2(nsolid),n3(nsolid),i,n1MAX,n2MAX,n3MAX
       double precision     :: rx(nsolid),ry(nsolid),rz(nsolid)
       double precision     :: a(3,3)
       double precision     :: a0,a02
       
       a02= a0/2.d0
       
         do i=1,nsolid
	        n1(i)=nint(rx(i)/a02)
	        n2(i)=nint(ry(i)/a02)
	        n3(i)=nint(rz(i)/a02)
         end do
	  n1MAX = maxval(n1)
	  n2MAX=  maxval(n2) 
	  n3MAX=  maxval(n3) 
	  
	  write(*,*) 'nsolid,a0',nsolid,a0
	  write(*,*), 'a',a
	  print*, 'n1MAX,n2MAX,n3MAX',n1MAX,n2MAX,n3MAX 
       
       return
       end 
!#####################################################################
       	  
       
       subroutine config_sia (nsolid,nsia,n1MAX,n2MAX,n3MAX,a0,&
                        ns1,ns2,ns3, &
                        rxloop2,ryloop2,rzloop2)
       
       implicit none
       integer  :: nsolid,nsia,n1MAX,n2MAX,n3MAX  
       integer  :: i,ns1(nsia),ns2(nsia),ns3(nsia)
       double precision  :: a0,a02
       double precision  :: rxloop2(nsia),ryloop2(nsia),rzloop2(nsia)     	
	a02=a0/2.0
          open (unit=33,file='no00_new.inp',status='unknown')
	  write(33,'(i7)') nsia
	    do i=1,nsia

	        ns1(i)=nint(rxloop2(i)/a02)
	        ns2(i)=nint(ryloop2(i)/a02)
	        ns3(i)=nint(rzloop2(i)/a02)
               write(33,'(3i7)') ns1(i),ns2(i),ns3(i) !,&
	    end do   
	     

	    do i=1,nsia
	      if ( (ns1(i) > n1MAX) .or.(ns2(i) > n2MAX) .or. &
	        (ns3(i) > n3MAX) ) then
		print*, 'Big problem in read_config_sia'
		stop
	      end if	    
	    end do 
	          
       return
       end 
!#####################################################################
       	  
       	  
       subroutine check_site (nsolid,nsia,n1,n2,n3,ns1,ns2,ns3,nsite)		         
 
       implicit none
       integer :: nsolid,nsia,itst,is,i
       integer :: n1(nsolid),n2(nsolid),n3(nsolid), &
                        nsite(nsia),ns1(nsia),ns2(nsia),ns3(nsia)
				
	       do i=1,4
!         write(*,*) n1(i),n2(i),n3(i)
       end do 

	
	
	
	
	itst=0
	do is=1,nsia
	   do i=1,nsolid
	       if (((n1(i)==ns1(is)).and.&
	            (n2(i)==ns2(is)).and.&
		     (n3(i)==ns3(is)))) then
		   itst =itst +1
		   nsite(is)= i
!         write(*,*) n1(i),n2(i),n3(i)
!         write(*,*) ns1(is),ns2(is),ns3(is)
	       end if
	    end do
	  end do     			
        if (itst.ne.nsia) then
	  print*, 'there is a big problem in input'
	  print*, 'Check check_site'
	  print*, itst, nsia
	  stop	   	
	end if	
				
	return
	end			       

!#####################################################################
       	  
       subroutine write_out (nall,nsolid,nsia,volc,nsite,rx,  &
                       ry,rz,rsx,rsy,rsz,a)
       
       implicit none
       integer  :: nall,nsolid,nsia,n,i,icard,in,icnt,is
       double precision     :: rx(nsolid),ry(nsolid),rz(nsolid),CELLDM
       integer  :: nsite(nsia)       
       double precision     :: rsx(2),rsy(2),rsz(2),x,y,z, volc,ascl        
       double precision,dimension(:,:) :: a(3,3)
  
       open(12,file='sia.xyz')
       open(13,file='unit_cell') 
         CELLDM=1.0
	 in=1
	 
	 ascl=MAXVAL(a*volc)

	open (15,file='formd_new.str',form='formatted',status='unknown')
	open (17,file='geom_new.data',form='formatted',status='unknown')
	open (18,file='formd_new_ndm_scl.str',form='formatted',status='unknown')
        
	write(15,'(i7)') nall
        write(15,'(3f10.5)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(15,'(3f10.5)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(15,'(3f10.5)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

	write(17,'(i7)') nall
	write(17,*) 1.0
        write(17,'(3f10.5)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(17,'(3f10.5)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(17,'(3f10.5)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

        write(18,'("1   1   1")') 
	write(18,'(3f10.5)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(18,'(3f10.5)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(18,'(3f10.5)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
	write(18,'(i8)') nall



10       format(I5,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
          write(13,*) 'number of atoms in the unit cell '
          write(13,*) nall 
          write(13,*) 'CELLDM= '
          write(13,*)  CELLDM
	  write(13,*) 'position of the atoms in the unit cell'



       write(12,'(i7)') nall
       write(12,'(a)') ' '
       icnt=0
       do i=1,nsia
        do n=1,2
	 icnt =icnt+1
         x = rx(nsite(i))*volc + rsx(n)*volc
         y = ry(nsite(i))*volc + rsy(n)*volc
         z = rz(nsite(i))*volc + rsz(n)*volc
       write(12,'("Fe  ",3f12.5)') x,y,z
       write(13,10) icnt,x,y,z,in,in,in
       write(15,'(a3,3(f10.5))') 'Fe ', x,y,z
       write(17,'(a3,3(f10.5))') '   ', x,y,z
       write(18,'(a3,3(f15.10),"   1")') '   ', x/ascl,y/ascl,z/ascl
	end do 
       end do
       
       
       do i=1,nsolid
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
       write(15,'(a3,3(f10.5))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
       write(17,'(a3,3(f10.5))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
       write(18,'(a3,3(f15.10),"   1")') '   ', rx(i)*volc/ascl,        &
                                        ry(i)*volc/ascl, rz(i)*volc/ascl
       end if
       end do
       
       write(12,'("  ")')
       
       write(12,'(3f10.5)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
       write(12,'(3f10.5)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
       write(12,'(3f10.5)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

       
       open(14,file='unit_vect',form='formatted') 
       write (14,*) ' number of unit vectors 0 1 2 ou 3' 
       write(14,*) 3
       write(14,*) 'unit vectors:'
       write(14,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
       write(14,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
       write(14,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
       

	


       close(12)
       close(13)
       close(14)
       


       return
       end


  subroutine VectProd(a1,a2,b)

    IMPLICIT NONE

    REAL(kind(0.d0)), dimension(1:3), intent(in) :: a1, a2
    REAL(kind(0.d0)), dimension(1:3) :: b

    b(1)=a1(2)*a2(3)-a1(3)*a2(2)
    b(2)=a1(3)*a2(1)-a1(1)*a2(3)
    b(3)=a1(1)*a2(2)-a1(2)*a2(1)

  return
  end

!#####################################################################
       
       
