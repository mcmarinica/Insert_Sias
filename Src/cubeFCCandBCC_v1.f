	  implicit none 
          integer :: istruct,n1,n2,n3,ntemp,natc, in,ic,itt,i,j,k
          real(8) :: a0,a02,CELLDM,Rx,Ry,Rz,x,y,z  
	  real(8) :: a(3,3), at(4,3)
	  
	  write(*,*) 'FCC=1/BCC=2'
	  read (*,*) istruct
	  write(*,*) 'Put the size of the cube Cu 3.615/ Fe 2.8553:'
	  read(*,*) a0
	  write(*,*) 'n1,n2,n3?'
	  read(*,*) n1,n2,n3
	  
	  a02=a0/2.d0
      open (10, file='formd.str', status='unknown',
     &      access='sequential', form='formatted')
	  
	  if (istruct.eq.1) then
	     natc=4
	     a(1,1)= a0
	     a(1,2)= 0.d0
	     a(1,3)= 0.d0
	     
	     a(2,1)= 0.d0
	     a(2,2)= a0
	     a(2,3)= 0.d0

	     a(3,1)= 0.d0
	     a(3,2)= 0.d0
	     a(3,3)= a0
	     
	     at(1,1)=0.d0
	     at(1,2)=0.d0
	     at(1,3)=0.d0
	     
	     at(2,1)=0
	     at(2,2)=a02
	     at(2,3)=a02

	     at(3,1)=a02
	     at(3,2)=a02
	     at(3,3)=0.d0
	     
	     at(4,1)=a02
	     at(4,2)=0.d0
	     at(4,3)=a02
	     
	     write(*,*) 'Writing FCC structure ...'
	    
	    else if (istruct.eq.2) then
	     
	     natc=2
	     a(1,1)= a0
	     a(1,2)= 0.d0
	     a(1,3)= 0.d0
	     
	     a(2,1)= 0.d0
	     a(2,2)= a0
	     a(2,3)= 0.d0

	     a(3,1)= 0.d0
	     a(3,2)= 0.d0
	     a(3,3)= a0
	     
	     at(1,1)=0.d0
	     at(1,2)=0.d0
	     at(1,3)=0.d0
	     
	     at(2,1)=a02
	     at(2,2)=a02
	     at(2,3)=a02
	     
	     write(*,*) 'Writing BCC structure ...'
	 end if    
	  
	  ntemp=n1*n2*n3*natc
	  
!writing the box
        write (*,'(3f10.5)') a(1,1),a(1,2), a(1,3)
        write (*,'(3f10.5)') a(2,1),a(2,2), a(2,3)
        write (*,'(3f10.5)') a(3,1),a(3,2), a(3,3)

	write (10,'(i8)') ntemp
        write (10,'(3f18.13)') 
     &	       a(1,1)*dble(n1), a(1,2)*dble(n1), a(1,3)*dble(n1)
        write (10,'(3f18.13)') 
     &	       a(2,1)*dble(n2), a(2,2)*dble(n2), a(2,3)*dble(n2)
        write (10,'(3f18.13)') 
     &	       a(3,1)*dble(n3), a(3,2)*dble(n3), a(3,3)*dble(n3)
	
       open(13,file='unit_cell',form='formatted') 
         CELLDM=1.0
	 in=1
10       format(I5,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
          write(13,*) 'number of atoms in the unit cell '
          write(13,*) ntemp 
          write(13,*) 'CELLDM= '
          write(13,*)  CELLDM
	  write(13,*) 'position of the atoms in the unit cell'

	        open(14,file='unit_vect',form='formatted') 
       write (14,*) ' number of unit vectors 0 1 2 ou 3' 
       write(14,*) 3
       write(14,*) 'unit vectors:'
       write(14,'(3f12.5)') 
     &          a(1,1)*dble(n1), a(1,2)*dble(n1), a(1,3)*dble(n1)
       write(14,'(3f12.5)') 
     &          a(2,1)*dble(n2), a(2,2)*dble(n2), a(2,3)*dble(n2)
       write(14,'(3f12.5)') 
     &          a(3,1)*dble(n3), a(3,2)*dble(n3), a(3,3)*dble(n3)
	
	ic=0
	itt=0
	do i=0,n1-1
	  do j=0,n2-1
	    do k=0,n3-1
	      ic=ic+1
	      Rx = a(1,1)*dble(i) + a(2,1)*dble(j) + a(3,1)*dble(k)
	      Ry = a(1,2)*dble(i) + a(2,2)*dble(j) + a(3,2)*dble(k)
	      Rz = a(1,3)*dble(i) + a(2,3)*dble(j) + a(3,3)*dble(k)
              do ic=1,natc
	        itt=itt+1
	        x = rx + at(ic,1)
	        y = ry + at(ic,2)
	        z = rz + at(ic,3)		
		
		write(10,'("Fe ",3(f18.13))')   x,y,z
	        write(11,'(i4,6(f10.5))') 1,x,y,z,0.d0,0.d0,0.d0
		write(13,10) itt,x,y,z,in,in,in

	      end do	
             end do
           end do
	 end do 
	 
	 
	 

	 
	end   	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      










	  
	  
	  
