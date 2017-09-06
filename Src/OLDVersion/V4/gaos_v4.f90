  program insert_c15_in_bcc
   use octahedre 
   implicit none
   real(8) :: alpha
   real(8) :: origin(3)
   integer :: i,j,iunit,iunitCu,iunitFe,type_octo

   call init_octa
      
   alpha=0.5d0
   iunit=21
   iunitCu=22
   iunitFe=23
   open(unit=iunit,name='temp.bs')
   open(unit=iunitCu,name='tempCu.bs')
   open(unit=iunitFe,name='tempFe.bs')

! astea baga in directia 110
   origin(1:3)=(/  1,  1,  1 /)  ! first type tetrahedron
   type_octo=1
   do i=1,4
      do j=1,3
        write(iunit,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+ origin(1:3)
        write(iunitFe,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+ origin(1:3)
      end do
        write(iunit,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
        write(iunitCu,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
   end do
   
   
   origin(1:3)=(/  2,  2,  2 /)  !rotated one
   type_octo=2
   do i=1,1
      do j=1,3
        write(iunit,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
        write(iunitFe,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
      end do
        write(iunit,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
        write(iunitCu,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
   end do

   do i=3,3
      do j=2,2
        write(iunit,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
        write(iunitFe,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
      end do
        write(iunit,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
        write(iunitCu,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
   end do


!d   origin(1:3)=(/  2,  0,  0 /)  !rotated
!d   type_octo=2
!d   do i=1,4
!d     do j=1,3
!d        write(iunit,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	  write(iunitFe,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	end do
!d	  write(iunit,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	  write(iunitCu,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
!d   end do
!d
!d   origin(1:3)=(/  0,  2,  0 /)  !rotated
!d   type_octo=2
!d   do i=1,4
!d	do j=1,3
!d	  write(iunit,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	  write(iunitFe,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	end do
!d	  write(iunit,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	  write(iunitCu,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
!d   end do
!d
!d! astea decaleaza si la fel pot sa continui in directia 110  baga in 
!d   origin(1:3)=(/  0,  0,  2 /)  !rotated one
!d   type_octo=2
!d   do i=1,4
!d	do j=1,3
!d	  write(iunit,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	  write(iunitFe,'("atom Fe   ",(3f8.3))') alpha*octo(type_octo)%gaos(1:3,j,i)+octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	end do
!d	  write(iunit,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
!d	  write(iunitCu,'("atom Cu   ",(3f8.3))') octo(type_octo)%atom(1:3,i)+origin(1:3)
!d   end do
!d
!d

   call print_suffix_bs(iunit)
   call print_suffix_bs(iunitFe)
   call print_suffix_bs(iunitCu)
   
   end program insert_c15_in_bcc
   
   
