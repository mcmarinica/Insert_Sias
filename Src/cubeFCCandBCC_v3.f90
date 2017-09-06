          !construct the cube by shells .......
	  ! you can and and substract shells to ab-initio
	  !
          implicit double precision (a-h,o-z)
	  
	  double precision a(3,3) 
	  real(8), dimension(:,:), allocatable :: at
	  write(*,*) 'FCC=1/BCC=2/TT=3'
	  read (*,*) istruct
	  write(*,*) 'Put the size of the cube Cu 3.61/ Fe 2.8553:'
	  read(*,*) a0
	  write(*,*) 'the outter shell n1,n2,n3?'
	  read(*,*) n1,n2,n3
	  write(*,*) 'the inner shell n1, n2,n3?'
	  read(*,*) nin1,nin2,nin3
	  
	  a02=a0/2.d0
   
	  if (istruct.eq.1) then
	     natc=4
             allocate(at(natc,3))
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
             allocate(at(natc,3))
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
	  

	    else if (istruct.eq.3) then
	     
	     natc=24
             allocate(at(natc,3))

	     a(1,1)= a0
	     a(1,2)= 0.d0
	     a(1,3)= 0.d0
	     
	     a(2,1)= 0.d0
	     a(2,2)= a0
	     a(2,3)= 0.d0

	     a(3,1)= 0.d0
	     a(3,2)= 0.d0
	     a(3,3)= a0
	     
             at( 1,:)=(/0.000d0,  0.000d0,  0.000d0/)
             at( 2,:)=(/0.250d0,  0.250d0,  0.250d0/)
             at( 3,:)=(/0.500d0,  0.500d0,  0.000d0/)
             at( 4,:)=(/0.500d0,  0.000d0,  0.500d0/)
             at( 5,:)=(/0.000d0,  0.500d0,  0.500d0/)
             at( 6,:)=(/0.750d0,  0.750d0,  0.250d0/)
             at( 7,:)=(/0.625d0,  0.125d0,  0.125d0/)
             at( 8,:)=(/0.750d0,  0.250d0,  0.750d0/)
             at( 9,:)=(/0.875d0,  0.375d0,  0.125d0/)
             at(10,:)=(/0.875d0,  0.125d0,  0.375d0/)
             at(11,:)=(/0.625d0,  0.375d0,  0.375d0/)
             at(12,:)=(/0.125d0,  0.625d0,  0.125d0/)
             at(13,:)=(/0.250d0,  0.750d0,  0.750d0/)
             at(14,:)=(/0.375d0,  0.875d0,  0.125d0/)
             at(15,:)=(/0.375d0,  0.625d0,  0.375d0/)
             at(16,:)=(/0.125d0,  0.875d0,  0.375d0/)
             at(17,:)=(/0.625d0,  0.625d0,  0.625d0/)
             at(18,:)=(/0.375d0,  0.375d0,  0.625d0/)
             at(19,:)=(/0.375d0,  0.125d0,  0.875d0/)
             at(20,:)=(/0.125d0,  0.375d0,  0.875d0/)
             at(21,:)=(/0.875d0,  0.875d0,  0.625d0/)
             at(22,:)=(/0.875d0,  0.625d0,  0.875d0/)
             at(23,:)=(/0.625d0,  0.875d0,  0.875d0/)
             at(24,:)=(/0.125d0,  0.125d0,  0.625d0/)

	     write(*,*) 'Writing TT structure ...'
             at(:,:)=a0*at(:,:)
	 end if    

      open (10, file='formd.str', status='unknown', &     
            access='sequential', form='formatted')
      open (15, file='formd.str_inner', status='unknown',  & 
           access='sequential', form='formatted')
      open (16, file='formd.gin', status='unknown',  & 
           access='sequential', form='formatted')


	
	  ntemp=n1*n2*n3*natc
	  
!writing the box
        write (*,'(3f10.5)') a(1,1:3)
        write (*,'(3f10.5)') a(2,1:3)
        write (*,'(3f10.5)') a(3,1:3)

	write (10,'(i6)') n1*n2*n3*natc
        write (10,'(3f18.13)') a(1,1:3)*dble(n1)
        write (10,'(3f18.13)') a(2,1:3)*dble(n2)
        write (10,'(3f18.13)') a(3,1:3)*dble(n3)

	write (15,'(i6)') nin1*nin2*nin3*natc
        write (15,'(3f18.13)')  a(1,1:3)*dble(nin1)
        write (15,'(3f18.13)')  a(2,1:3)*dble(nin2)
        write (15,'(3f18.13)')  a(3,1:3)*dble(nin3)

        write (16,*) "1 1 1"
        write (16,'(3f18.13)') a(1,1:3)*dble(n1)
        write (16,'(3f18.13)') a(2,1:3)*dble(n2)
        write (16,'(3f18.13)') a(3,1:3)*dble(n3)

	write (16,'(i6)') n1*n2*n3*natc

	

	ic=0
	itt=0
	do i=0,nin1-1
	  do j=0,nin2-1
	    do k=0,nin3-1
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
		write(15,'("Fe ",3(f18.13))')   x,y,z
	        write(11,'(i4,6(f10.5))')1,x,y,z,0,0,0

	        write(16,'(3(f25.15), "  1")')x/(a(1,1)*dble(n1)),y/(a(2,2)*dble(n2)) & 

                                       ,z/(a(3,3)*dble(n3))
	      end do	
             end do
           end do
	 end do 
	 
	do i=0,n1-1
	  do j=0,n2-1
	    do k=0,n3-1
	     if ( (i.ge.nin1).or.(j.ge.nin2).or.(k.ge.nin3) ) then 
	      ic=ic+1
	      Rx = a(1,1)*dble(i) + a(2,1)*dble(j) + a(3,1)*dble(k)
	      Ry = a(1,2)*dble(i) + a(2,2)*dble(j) + a(3,2)*dble(k)
	      Rz = a(1,3)*dble(i) + a(2,3)*dble(j) + a(3,3)*dble(k)
              do ic=1,natc
	        itt=itt+1
	        x = Rx + at(ic,1)
	        y = Ry + at(ic,2)
	        z = Rz + at(ic,3)		
		
	
	write(10,'("Fe ",3(f18.13))')   x,y,z
	        write(11,'(i4,6(f10.5))')1,x,y,z,0,0,0
	        write(16,'(3(f25.15), "  1")')x/(a(1,1)*dble(n1)),y/(a(2,2)*dble(n2)) & 
                                       ,z/(a(3,3)*dble(n3))

	      end do	
	    end if  
             end do
           end do
	 end do 
	 
	 

	 
	end   	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      
	      










	  
	  
	  
