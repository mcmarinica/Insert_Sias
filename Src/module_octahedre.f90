module octahedre

   implicit none
   real(8) :: pi=2.d0*acos(0.d0)
   real(8) :: deg2rad=2.d0*acos(0.d0)/180.d0
   type octhaedre   
   real(8) :: center(3)    ! center of the octhaedre
   real(8) :: link(3,4)    ! the first direction which can link an other
                           ! tetrahedron ( can be the the center of an 
			   ! octa from the second type)
   real(8) :: atom(3,4)          ! vacancy in the corner
   real(8) :: removeatom(3,6)    ! vacancy up-down side 
   real(8) :: gaos(3,3,4)  ! first index x,y,z of the one of the three gao atoms
                           ! aroud the atom 1 to 4
   end type octhaedre

   type(octhaedre) , dimension(4) :: octo

contains
 
  subroutine init_octa
  integer i,j
   real(8) :: vec_temp_in(3),vec_temp_out(3)
!........construction of the octo's............
   !the first octo is based in the first cube 111 center
   
   !the coordinates of the  tetrahedron
   octo(1)%atom(1:3,1)=(/ -1, -1, -1 /)  !000   
   octo(1)%atom(1:3,2)=(/  1,  1, -1 /)  !220  
   octo(1)%atom(1:3,3)=(/  1, -1,  1 /)  !202 
   octo(1)%atom(1:3,4)=(/ -1,  1,  1 /)  !022
   
   !the other six atoms which must be removed:
   octo(1)%removeatom(1:3,1)=(/  0,  0,  -2 /)
   octo(1)%removeatom(1:3,2)=(/  0,  2,   0 /)
   octo(1)%removeatom(1:3,3)=(/  0,  0,   2 /)
   octo(1)%removeatom(1:3,4)=(/  0, -2,   0 /)
   octo(1)%removeatom(1:3,5)=(/ -2,  0,   0 /)
   octo(1)%removeatom(1:3,6)=(/  2,  0,   0 /)
   
   !so be wise budy ...you remove 10 (4 in the corner of the tetrahedron 
   ! and 6 up-down side) atomes and you  put 12 (grouped in the 4 gaos). 
   !Finally you will have an di-interstitials  
   
   !the 4 directions which can be linked (the atom coordinates ca be obtained
   ! using  octo(1)%link(1:3,i) +  octo(1)%center(1:3), i=1,4
   octo(1)%link(1:3,1)=  (/  1,   1 ,   1 /)
   octo(1)%link(1:3,2)=  (/ -1,  -1 ,   1 /)
   octo(1)%link(1:3,3)=  (/  1,  -1 ,  -1 /)
   octo(1)%link(1:3,4)=  (/ -1,   1 ,  -1 /)

   octo(1)%gaos(1:3,1,1)=(/  1 ,  1 ,  -1 /)
   octo(1)%gaos(1:3,2,1)=(/  1 , -1 ,   1 /)
   octo(1)%gaos(1:3,3,1)=(/ -1 ,  1 ,   1 /) 

   octo(1)%gaos(1:3,1,2)=(/ -1 , -1 ,  -1 /)
   octo(1)%gaos(1:3,2,2)=(/ -1 ,  1 ,   1 /)
   octo(1)%gaos(1:3,3,2)=(/  1 , -1 ,   1 /) 

   octo(1)%gaos(1:3,1,3)=(/ -1 ,  1 ,   1 /)
   octo(1)%gaos(1:3,2,3)=(/ -1 , -1 ,  -1 /)
   octo(1)%gaos(1:3,3,3)=(/  1 ,  1 ,  -1 /) 

   octo(1)%gaos(1:3,1,4)=(/  1 , -1 ,   1 /)
   octo(1)%gaos(1:3,2,4)=(/  1 ,  1 ,  -1 /)
   octo(1)%gaos(1:3,3,4)=(/ -1 , -1 ,  -1 /) 

   octo(1)%center(1:3)=(/  1,  1,  1 /)  
   
   !The second cube has the centere in 222 and is rotated 90 counterclockwise and 222 based
   octo(2)%center(1:3)=(/  2,  2,  2 /)
   do i=1,4
      vec_temp_in(1:3)=octo(1)%atom(1:3,i)
      !debug write(*,*) 'octo',vec_temp_in
       call rotation_vector_z(90*deg2rad,vec_temp_in,vec_temp_out)
      !debug write(*,*) 'octo',vec_temp_out
      octo(2)%atom(1:3,i)=vec_temp_out(1:3)
      
      do j=1,3
       vec_temp_in(1:3)=octo(1)%gaos(1:3,j,i)
       call rotation_vector_z(90*deg2rad,vec_temp_in,vec_temp_out)
       octo(2)%gaos(1:3,j,i)=vec_temp_out(1:3)       
      end do
                  
   end do      

      do j=1,4
      vec_temp_in(1:3)=octo(1)%link(1:3,j)
       call rotation_vector_z(90*deg2rad,vec_temp_in,vec_temp_out)
       octo(2)%link(1:3,j)=vec_temp_out(1:3)       
      end do
   
      do j=1,6
      vec_temp_in(1:3)=octo(1)%removeatom(1:3,j)
       call rotation_vector_z(90*deg2rad,vec_temp_in,vec_temp_out)
       octo(2)%removeatom(1:3,j)=vec_temp_out(1:3)       
      end do

  return
  end  subroutine init_octa



 end module octahedre

module icosaedre

   implicit none
   integer, dimension(:), allocatable :: type_icos
   real(8) :: pi=2.d0*acos(0.d0)
   real(8) :: deg2rad=2.d0*acos(0.d0)/180.d0
   real(8), dimension(3)  :: t_dir

   type icosah 
   real(8) :: center(3)    ! center of the octhaedre
   real(8) :: link(3,4)    ! the first direction which can link an other
                           ! tetrahedron ( can be the the center of an 
			   ! octa from the second type)
   real(8) :: sia(3,12)          ! vacancy in the corner
   real(8) :: add(3,1)          ! vacancy in the corner
   real(8) :: removeatom(3,6)    ! vacancy up-down side 
                           ! aroud the atom 1 to 4
   end type icosah

   type(icosah) , dimension(4) ::  icos

contains
 
  subroutine init_icos
  integer i,j
   real(8) :: vec_temp_in(3),vec_temp_out(3)
!........construction of the octo's............
   !the first octo is based in the first cube 111 center
   
   
   !the other 8 atoms which must be removed (each one one face):

   icos(1)%removeatom(1:3,1)=(/ -1.d0,  0.d0,   0.d0 /) !011
   icos(1)%removeatom(1:3,2)=(/  1.d0,  0.d0,   0.d0 /) !211
   icos(1)%removeatom(1:3,3)=(/  0.d0,  0.d0,  -1.d0 /) !110
   icos(1)%removeatom(1:3,4)=(/  0.d0,  0.d0,   1.d0 /) !112
   icos(1)%removeatom(1:3,5)=(/  0.d0, -1.d0,   0.d0 /) !101
   icos(1)%removeatom(1:3,6)=(/  0.d0,  1.d0,   0.d0 /) !121
   

   !atom in the center 
   icos(1)%add(1:3,1)=(/  0.0, 0.0,   0.0 /)  !111

   !the other 12 atoms in SIAs position:
   icos(1)%sia(1:3, 1)=(/ -1.0d0,  0.0d0,   0.5d0 /) !011
   icos(1)%sia(1:3, 2)=(/ -1.0d0,  0.0d0,  -0.5d0 /) !011

   icos(1)%sia(1:3, 3)=(/  1.0d0,  0.0d0,   0.5d0 /) !211
   icos(1)%sia(1:3, 4)=(/  1.0d0,  0.0d0,  -0.5d0 /) !211

   icos(1)%sia(1:3, 5)=(/  0.0d0,  0.5d0,  -1.0d0 /) !110
   icos(1)%sia(1:3, 6)=(/  0.0d0, -0.5d0,  -1.0d0 /) !110

   icos(1)%sia(1:3, 7)=(/  0.0d0,  0.5d0,   1.0d0 /) !112
   icos(1)%sia(1:3, 8)=(/  0.0d0, -0.5d0,   1.0d0 /) !112

   icos(1)%sia(1:3, 9)=(/  0.5d0, -1.0d0,   0.0d0 /) !101
   icos(1)%sia(1:3,10)=(/ -0.5d0, -1.0d0,   0.0d0 /) !101


   icos(1)%sia(1:3,11)=(/  0.5d0,  1.0d0,   0.0d0 /) !121
   icos(1)%sia(1:3,12)=(/ -0.5d0,  1.0d0,   0.0d0 /) !121



   !atom in the center 
   icos(2)%add(1:3,1)=(/  0.0, 0.0,   0.0 /)  !111
   t_dir=(/ 1.d0/dsqrt(3.d0),  1.d0/dsqrt(3.d0), 1.d0/dsqrt(3.d0) /)
   do i=1,6
      vec_temp_in(1:3)=icos(1)%removeatom(1:3,i)
      !debug write(*,*) 'octo',vec_temp_in
       call rotation_vector_t(60*deg2rad,t_dir, vec_temp_in,vec_temp_out)
       write(*,*) 'icos',vec_temp_in, 't', vec_temp_out
      icos(2)%removeatom(1:3,i)=vec_temp_out(1:3)
   end do      

   do j=1,12
     vec_temp_in(1:3)=icos(1)%sia(1:3,j)
     call rotation_vector_t(60*deg2rad,t_dir, vec_temp_in,vec_temp_out)
     icos(2)%sia(1:3,j)=vec_temp_out(1:3)       
   end do

  return
  end  subroutine init_icos



 end module icosaedre




