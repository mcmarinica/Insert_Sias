


!#####################################################################
	 	
	
       subroutine read_ntemp (ntemp)
       integer :: ntemp
        
	open (10,file='formd.str',form='formatted',status='unknown')
        read (10,'(i6)') ntemp
       
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
       	  
       subroutine read_input(nall,ntemp,a,rx,ry,rz)      
       integer              :: nall, ntemp
       double precision     :: rx(nall),ry(nall),rz(nall)
       double precision,dimension(:,:) :: a(3,3)
       character*3 :: sym
       
        read(10,'(3f10.5)')a(1,1), a(1,2), a(1,3)
        read(10,'(3f10.5)')a(2,1), a(2,2), a(2,3)
        read(10,'(3f10.5)')a(3,1), a(3,2), a(3,3)
	
	do i=1,ntemp
		read(10,'(a3,3(f10.5))') sym, rx(i),ry(i),rz(i)
        end do 		
     
		
       
       return
       end 
!#####################################################################
