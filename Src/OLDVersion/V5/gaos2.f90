   implicit none
   
   real(8) :: pi=2.d0*acos(0.d0)
   real(8) :: deg2rad=2.d0*acos(0.d0)/180.d0
   real(8) :: vec_temp_in(3),vec_temp_out(3)
   real(8) :: matrix_in(3,3),matrix_out(3,3),alpha
   integer :: i,i1,i2,i3,j
   type octhaedre   
   real(8) :: center(3)
   real(8) :: atom(3,4)
   end type octhaedre

   real(8)  :: mat_gaos(3,3)
   type(octhaedre) , dimension(2) :: octo
   
!........construction of the octo's............
   !the first octo is based in the first cube 111 center
   octo(1)%atom(1:3,1)=(/ -1, -1, -1 /)   
   octo(1)%atom(1:3,2)=(/  1,  1, -1 /)   
   octo(1)%atom(1:3,3)=(/  1, -1,  1 /)   
   octo(1)%atom(1:3,4)=(/ -1,  1,  1 /)

   octo(1)%center(1:3)=(/  1,  1,  1 /)  
    
   !The second cube is rotated 90 counterclockwise and 222 based
   octo(2)%center(1:3)=(/  2,  2,  2 /)
   do i=1,4
      vec_temp_in(1:3)=octo(1)%atom(1:3,i)
       call rotation_vector_z(90*deg2rad,vec_temp_in,vec_temp_out)
      octo(2)%atom(1:3,i)=vec_temp_out(1:3)
   end do   
   
   write(*,*) 'First octo'
   do i=1,4
    write(*,'(3f4.1)') octo(1)%atom(1:3,i)+octo(1)%center(1:3)
   end do
   

   write(*,*) 'Second octo'
   do i=1,4
    write(*,'(3f4.1)') octo(2)%atom(1:3,i)+octo(2)%center(1:3)
   end do
!......end construction of the gao's....................   
   
   
   mat_gaos(1:3,1)=(/  1 ,  1 ,  -1 /)
   mat_gaos(1:3,2)=(/  1 , -1 ,   1 /)
   mat_gaos(1:3,3)=(/ -1 ,  1 ,   1 /) 
   
   alpha=0.25d0
   

! construction of the first tetrahedron
   do i=1,4
   i1= int(octo(1)%atom(1,i) + octo(1)%center(1))
   i2= int(octo(1)%atom(2,i) + octo(1)%center(2))
   i3= int(octo(1)%atom(3,i) + octo(1)%center(3))
   
   matrix_in(:,:)=mat_gaos(:,:)
!   do j=1,3
!   write(*,'("g",3f6.2)') mat_gaos(1:3,j)
!   end do
   write(*,*) i1,i2,i3
   call reflexion(1,i1,matrix_in,matrix_out)

!   write(*,'("1",3f6.2)') alpha*matrix_out(1:3,1) 
!   write(*,'("1",3f6.2)') alpha*matrix_out(1:3,2) 
!   write(*,'("1",3f6.2)') alpha*matrix_out(1:3,3) 

   call reflexion(2,i2,matrix_out,matrix_in)
   call reflexion(3,i3,matrix_in,matrix_out)
   
!   write(*,*) '------------------------'
   write(*,'(3f6.2)') alpha*matrix_out(1:3,1) 
   write(*,'(3f6.2)') alpha*matrix_out(1:3,2) 
   write(*,'(3f6.2)') alpha*matrix_out(1:3,3) 
   
   end do
    
   
   end
   
   
   
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
  real(8)  :: pi=2.d0*acos(0.d0)
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
  real(8)  :: pi=2.d0*acos(0.d0)
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
