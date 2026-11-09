!=======================================================================
! rotations.f90
!
! Rotation of a vector, used to build the second orientation of the C15
! and icosahedral units.
!
!   rotation_vector_z  rotation about z
!   rotation_vector_t  rotation about an arbitrary axis
!=======================================================================

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

!=======================================================================

  subroutine rotation_vector_t(teta, tdir, vector_in,vector_out)
  !Rotation along the direction tdir with the angle teta
  implicit none
  real(8), intent(in) :: teta
  real(8), dimension(3), intent(in) :: tdir
  real(8), dimension(3), intent(in)  :: vector_in
  real(8), dimension(3), intent(out) :: vector_out
  real(8) :: atemp, rcos, rsin
  integer :: i,k
  real(8), dimension(3,3) :: mat_rotation
  rcos =  dcos(teta)
  rsin =  dsin(teta)

  mat_rotation(1,1:3)=(/  tdir(1)**2*(1.d0-rcos)     +        rcos, &
                          tdir(1)*tdir(2)*(1.d0-rcos)-tdir(3)*rsin, &
                          tdir(1)*tdir(3)*(1.d0-rcos)+tdir(2)*rsin  /)
  mat_rotation(2,1:3)=(/  tdir(1)*tdir(2)*(1.d0-rcos)+tdir(3)*rsin, &
                          tdir(2)**2*(1.d0-rcos)     +        rcos, &
                          tdir(2)*tdir(3)*(1.d0-rcos)-tdir(1)*rsin  /)
  mat_rotation(3,1:3)=(/  tdir(1)*tdir(3)*(1.d0-rcos)-tdir(2)*rsin, &
                          tdir(2)*tdir(3)*(1.d0-rcos)+tdir(1)*rsin, &
                          tdir(3)**2*(1.d0-rcos)     +        rcos  /)


    do i=1,3
	atemp=0
        do k=1,3
	  atemp=atemp+mat_rotation(i,k)*vector_in(k)
	end do 
      vector_out(i)=atemp
     end do
  return
  end

!=======================================================================
