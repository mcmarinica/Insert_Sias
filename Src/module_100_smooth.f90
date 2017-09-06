 module the_100_smooth

implicit none
 integer :: pbc_100
 real(8) :: rc100x(-2:2),rc100y(-2:2),rc100z(-2:2)
 integer, allocatable, dimension(:,:) :: nsite_100,ns1_100,ns2_100,ns3_100

 contains 


 subroutine allocate_the_100 (nsia)
 implicit none
 integer, intent(in) :: nsia
 
 allocate (nsite_100(-2:2,nsia),ns1_100(-2:2,nsia),ns2_100(-2:2,nsia),ns3_100(-2:2,nsia))

 return

end subroutine allocate_the_100

 subroutine fill_the_100 (a0,as100,rsx,rsy,rsz)

 implicit none
 real(8), intent(in) :: a0,as100
 real(8), dimension(2), intent (in) :: rsx,rsy,rsz
 real(8) :: ann,dx,dx1,dx2,x_low_limit,x_high_limit
! Fixing the limits. We start from the assumption the rsx(1)=rsx(2)=rsz(3)=d/2
!The upper limits of x is done by:

   

   ann=a0
   if ((1.5d0*dabs(as100) - ann)<=0) then
    write(*,*) 'as100 too small compared to a0. Increase as100'
    write(*,*) 'stop in <fill_the_100>'
    stop 
   end if

   if (as100>=(ann*dble(pbc_100)/(dble(pbc_100)+1.d0))) then
     write(*,*) 'Reduce as100 to be smaller than ...:',ann*dble(pbc_100)/(dble(pbc_100)+1.d0)
     write(*,*) 'Stop in <fill_the_100>'
     stop
   end if
   x_low_limit=1.5d0*as100 - ann
   x_high_limit=ann*dble(pbc_100)/(dble(pbc_100)+1)-ann+0.5d0*as100
   dx= 0.5d0*(x_low_limit + x_high_limit)
   write(*,*) 'HELP', dx, x_low_limit, x_high_limit

   
    if (dx < as100/3.d0) then
       dx1=dx
       dx2=0.5d0*(5d0*dx/2.d0-as100/2.d0)
    end if       

   
    if (dx >= as100/3.d0) then
      if (as100/3.d0<=x_high_limit) then
       dx= 0.5d0*(dx + x_low_limit)
       dx1=dx
       dx2=0.5d0*(5d0*dx/2.d0-as100/2.d0)
      end if 
      if (as100/3.d0 > x_high_limit ) then 
        write(*,*) 'No solution for this value of as111', as100 
        write(*,*) 'Choose something in between ...', ann/1.5d0 ,' and ', ann*dble(pbc_100)/(dble(pbc_100)+1.d0)
        write(*,*) 'but no bigger than ...', x_high_limit/3.d0
        write(*,*) 'but no smaller than ...', x_low_limit/3.d0
        stop
      end if  
    end if  


 
    rc100x(-2)=0.d0
    rc100y(-2)=0.d0
    rc100z(-2)=-dabs(dx2)

    rc100x(-1)=0.d0
    rc100y(-1)=0.d0
    rc100z(-1)=-dabs(dx1)

    rc100x(0)=0.d0
    rc100y(0)=0.d0
    rc100z(0)=0.d0

    rc100x(1)=-rc100x(-1)
    rc100y(1)=-rc100y(-1)
    rc100z(1)=-rc100z(-1)

    rc100x(2)=-rc100x(-2)
    rc100y(2)=-rc100y(-2)
    rc100z(2)=-rc100z(-2)
 
 return
end subroutine fill_the_100



 subroutine build_the_nsite_the_100(nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
 implicit none
  integer, intent(in) :: nsia,n1MAX,n2MAX,n3MAX
  integer, intent (in) :: ns1(nsia),ns2(nsia),ns3(nsia)
  integer  :: is,jj,n1pbc,n2pbc,n3pbc

  n1pbc=n1MAX+1
  n2pbc=n2MAX+1
  n3pbc=n3MAX+1


 do is=1,nsia
  do jj=-2,2,1
  ns1_100(jj,is)=ns1(is)
  ns2_100(jj,is)=ns2(is)
  ns3_100(jj,is)=ns3(is)+2*jj
  end do
 end do



 do is=1,nsia
  do jj=-2,2,1
   if (ns1_100(jj,is)< 0) ns1_100(jj,is)= ns1_100(jj,is)+n1pbc
   if (ns2_100(jj,is)< 0) ns2_100(jj,is)= ns2_100(jj,is)+n2pbc
   if (ns3_100(jj,is)< 0) ns3_100(jj,is)= ns3_100(jj,is)+n3pbc

   if (ns1_100(jj,is)>= n1pbc ) ns1_100(jj,is)= ns1_100(jj,is)-n1pbc
   if (ns2_100(jj,is)>= n2pbc ) ns2_100(jj,is)= ns2_100(jj,is)-n2pbc
   if (ns3_100(jj,is)>= n3pbc ) ns3_100(jj,is)= ns3_100(jj,is)-n3pbc
 end do
end do

 return

end subroutine  build_the_nsite_the_100



subroutine distord_lattice_100 (ntemp,nsia,rx,ry,rz)
implicit none
integer , intent(in) :: nsia,ntemp
real(8), intent(inout) :: rx(ntemp),ry(ntemp),rz(ntemp)
integer :: jj,is,i

 do jj=-2,2,1
  do is=1,nsia
     i=nsite_100(jj,is)
     rx(i) = rx(i)+rc100x(jj)
     ry(i) = ry(i)+rc100y(jj)
     rz(i) = rz(i)+rc100z(jj)
     write(*,'(i5,3f12.5)') i, rc100x(jj),rc100y(jj),rc100z(jj)
   end do
  end do

 return

end subroutine distord_lattice_100


end module the_100_smooth


