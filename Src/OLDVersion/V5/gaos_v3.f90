  program insert_c15_in_bcc
   use octahedre 
   implicit none
   real(8) :: alpha
   real(8) :: origin(3)
   integer :: i,j,iunit

   call init_octa
      
   alpha=0.5d0
   iunit=21
   open(unit=iunit,name='temp.bs')

! astea baga in directia 110
   origin(1:3)=(/  1,  1,  1 /)
   do i=1,4
      do j=1,3
        write(21,'("atom Fe   ",(3f8.3))') alpha*octo(1)%gaos(1:3,j,i)+octo(1)%atom(1:3,i)+ origin(1:3)
      end do
        write(21,'("atom Cu   ",(3f8.3))') octo(1)%atom(1:3,i)+origin(1:3)
   end do
   
   
   origin(1:3)=(/  2,  2,  2 /)  !rotated one
   do i=1,4
      do j=1,3
        write(21,'("atom Fe   ",(3f8.3))') alpha*octo(2)%gaos(1:3,j,i)+octo(2)%atom(1:3,i)+origin(1:3)
      end do
        write(21,'("atom Cu   ",(3f8.3))') octo(2)%atom(1:3,i)+origin(1:3)
   end do


   origin(1:3)=(/  3,  3,  1 /)
   do i=1,4
      do j=1,3
        write(21,'("atom Fe   ",(3f8.3))') alpha*octo(1)%gaos(1:3,j,i)+octo(1)%atom(1:3,i)+origin(1:3)
      end do
        write(21,'("atom Cu   ",(3f8.3))') octo(1)%atom(1:3,i)+origin(1:3)
   end do

   origin(1:3)=(/  4,  4,  2 /)  !rotated
   do i=1,4
      do j=1,3
        write(21,'("atom Fe   ",(3f8.3))') alpha*octo(2)%gaos(1:3,j,i)+octo(2)%atom(1:3,i)+origin(1:3)
      end do
        write(21,'("atom Cu   ",(3f8.3))') octo(2)%atom(1:3,i)+origin(1:3)
   end do

! astea decaleaza si la fel pot sa continui in directia 110  baga in 
   origin(1:3)=(/  2,  0,  0 /)  !rotated one
   do i=1,4
      do j=1,3
        write(21,'("atom Fe   ",(3f8.3))') alpha*octo(2)%gaos(1:3,j,i)+octo(2)%atom(1:3,i)+origin(1:3)
      end do
        write(21,'("atom Cu   ",(3f8.3))') octo(2)%atom(1:3,i)+origin(1:3)
   end do



   call print_suffix_bs(iunit)
   
   end program insert_c15_in_bcc
   
   
   subroutine reflexion(type_r,n_reflexion,matrix_in,matrix_out)
   !
   implicit none
   integer,intent(in) :: type_r,n_reflexion
   real(8),intent(in)    ::matrix_in(3,3) 
   real(8),intent(out)   ::matrix_out(3,3)
   integer :: i,j
   real(8) :: factor
   
    if (type_r.gt.3) then
     write(*,*) 'ERROR .... no defined reflexion'
     write(*,*) 'type_r...........',type_r
    end if 
     do j=1,3
        factor=1.d0 
        if (j.eq.type_r) factor=(-1.d0)**n_reflexion
      do i=1,3
       matrix_out(j,i)=factor*matrix_in(j,i)     
      end do
    end do  
  return
  end
  
  subroutine rotation_matrix_z(teta,matrix_in,matrix_out)
  
  implicit none
  real(8) :: teta,atemp
  integer  :: i,j,k
  real(8), dimension(3,3) :: matrix_in,matrix_out,mat_rotation
  
  mat_rotation(1,1:3)=(/  dcos(teta),  -dsin(teta), 0.d0  /)
  mat_rotation(2,1:3)=(/  dsin(teta),   dcos(teta), 0.d0  /)
  mat_rotation(3,1:3)=(/        0.d0,         0.d0, 1.d0  /)
    do i=1,3
     do j=1,3 
	atemp=0
        do k=1,3
	  atemp=atemp+mat_rotation(i,k)*matrix_in(k,j)
	end do 
      matrix_out(i,j)=atemp
     end do
    end do
  return
  end

  subroutine rotation_vector_z(teta,vector_in,vector_out)
  
  implicit none
  real(8) :: teta,atemp
  integer :: i,k
  real(8), dimension(3) :: vector_in,vector_out
  real(8), dimension(3,3) :: mat_rotation
  
  mat_rotation(1,1:3)=(/  dcos(teta),  -dsin(teta), 0.d0  /)
  mat_rotation(2,1:3)=(/  dsin(teta),   dcos(teta), 0.d0  /)
  mat_rotation(3,1:3)=(/        0.d0,         0.d0, 1.d0  /)


    do i=1,3
	atemp=0
        do k=1,3
	  atemp=atemp+mat_rotation(i,k)*vector_in(k)
	end do 
      vector_out(i)=atemp
     end do
  return
  end
  
  subroutine print_suffix_bs(iunit)
  
    integer,intent(in) :: iunit
    
    
    write(iunit,*) 'spec   Fe    0.1000    0.6000'                                     
    write(iunit,*) 'spec   Fe2    0.1000    0.6000 					       '
    write(iunit,*) 'spec  Cu 0.1 0.1 0.5								   '
    write(iunit,*) '									   '
    write(iunit,*) 'bonds	Fe Fe	 0.0000    1.8000    0.0200    0.2000			   '
    write(iunit,*) 'bonds	Cu Cu	 2.1000    2.9000    0.0100    0.1 0.1 0.5		   '
    write(iunit,*) 'bonds	Fe Fe2    0.0000    0.8000    0.0200	0.1 0.1 0.9		   '
    write(iunit,*) '									   '
    write(iunit,*) 'tmat  1.0000  0.0000  0.0000  0.0000  0.0000  1.0000  0.0000 -1.0000  0.0000  '
    write(iunit,*) 'dist	 15.0000							   '
    write(iunit,*) 'inc	10.0000 							   '
    write(iunit,*) 'scale	  90.0000							   '
    write(iunit,*) 'rfac	  1.0000							   '
    write(iunit,*) 'bfac	  1.4000							   '
    write(iunit,*) 'pos	 0.0000      0.0000						   '
    write(iunit,*) 'light	   0.0000      0.0000	   0.0000				   '
    write(iunit,*) 'switches 1 0 1 0 0 1 1 0 0                                                    '
  
  
  return
  end
