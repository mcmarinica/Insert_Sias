subroutine define_center (p,n1MAX,surface,plane_corner)
implicit none
integer, intent(in)                 :: p,n1MAX
real(8), intent(out)                :: surface
real(8), dimension(3,4), intent(out) :: plane_corner
integer :: n2


! We cut in the plane (0,0,1). Here are the projection on this plane.
n2=(n1MAX-1)/2
plane_corner(:,1)= dble((/               n2,          0,     0 /))
plane_corner(:,2)= dble((/        n1MAX-1-1,       n2-1,     0 /))
plane_corner(:,3)= dble((/  n1MAX-1-1-(p-1),   n2-1+p-1,     0 /))
plane_corner(:,4)= dble((/         n2-(p-1),        p-1,     0 /))


surface= dsqrt(SUM((plane_corner(:,2)-plane_corner(:,1))**2))* &
         dsqrt(SUM((plane_corner(:,4)-plane_corner(:,1))**2))

!debug write(*,*) '......', dsqrt(SUM((plane_corner(:,2)-plane_corner(:,1))**2))
!debug write(*,*) '......', dsqrt(SUM((plane_corner(:,4)-plane_corner(:,1))**2))



return

end subroutine define_center
