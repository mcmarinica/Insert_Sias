 module the_111_smooth

implicit none
 integer :: pbc_111
 real(8) :: rcx(-2:2),rcy(-2:2),rcz(-2:2)
 !real(8), allocatable,dimension(:) :: deb_dx(:),deb_dy(:),deb_dz(:)
 !real(8), allocatable,dimension(:) :: fin_dx(:),fin_dy(:),fin_dz(:)
 integer, allocatable, dimension(:,:) :: nsite_111,ns1_111,ns2_111,ns3_111
 integer, allocatable, dimension(:,:) :: deb_nsite_111,deb_ns1_111,deb_ns2_111,deb_ns3_111
 integer, allocatable, dimension(:,:) :: fin_nsite_111,fin_ns1_111,fin_ns2_111,fin_ns3_111

 contains 


 subroutine allocate_the_111 (nsia)
 implicit none
 integer, intent(in) :: nsia
 
 allocate (nsite_111(-2:2,nsia),ns1_111(-2:2,nsia),ns2_111(-2:2,nsia),ns3_111(-2:2,nsia))

 return

end subroutine allocate_the_111

 subroutine allocate_the_111_migration (ntemp,nsia)
 implicit none
 integer, intent(in) :: nsia,ntemp
 
 allocate (deb_nsite_111(-2:2,nsia),deb_ns1_111(-2:2,nsia),deb_ns2_111(-2:2,nsia),deb_ns3_111(-2:2,nsia))
 allocate (fin_nsite_111(-2:2,nsia),fin_ns1_111(-2:2,nsia),fin_ns2_111(-2:2,nsia),fin_ns3_111(-2:2,nsia))
 !allocate (deb_dx(ntemp),deb_dy(ntemp),deb_dz(ntemp),fin_dx(ntemp),fin_dy(ntemp),fin_dz(ntemp))
 return

end subroutine allocate_the_111_migration


 subroutine fill_the_111 (a0,as111,rsx,rsy,rsz)

 implicit none
 real(8), intent(in) :: a0,as111
 real(8), dimension(2), intent (in) :: rsx,rsy,rsz
 real(8) :: ann,dx,dx1,dx2,x_low_limit,x_high_limit
! Fixing the limits. We start from the assumption the rsx(1)=rsx(2)=rsz(3)=d/2
!The upper limits of x is done by:

   

   ann=a0*dsqrt(3.d0)/2.d0
   if ((1.5d0*dabs(as111) - ann)<=0) then
    write(*,*) 'as111 too small compared to a0. Increase as111'
    write(*,*) 'stop in <fill_the_111>'
    stop 
   end if

   if (as111>=(ann*dble(pbc_111)/(dble(pbc_111)+1.d0))) then
     write(*,*) 'Reduce as111 to be smaller than ...:',ann*dble(pbc_111)/(dble(pbc_111)+1.d0)
     write(*,*) 'Stop in <fill_the_111>'
     stop
   end if
   x_low_limit=1.5d0*as111 - ann
   x_high_limit=ann*dble(pbc_111)/(dble(pbc_111)+1)-ann+0.5d0*as111
   dx= 0.5d0*(x_low_limit + x_high_limit)
   write(*,*) 'HELP', dx, x_low_limit, x_high_limit

   
    if (dx < as111/3.d0) then
       dx1=dx/dsqrt(3.d0)
       dx2=0.5d0*(5d0*dx/2.d0-as111/2.d0)/dsqrt(3.d0)
    end if       

   
    if (dx >= as111/3.d0) then
      if (as111/3.d0<=x_high_limit) then
       dx= 0.5d0*(dx + x_low_limit)
       dx1=dx/dsqrt(3.d0)
       dx2=0.5d0*(5d0*dx/2.d0-as111/2.d0)/dsqrt(3.d0)
      end if 
      if (as111/3.d0 > x_high_limit ) then 
        write(*,*) 'No solution for this value of as111', as111 
        write(*,*) 'Choose something in between ...', ann/1.5d0 ,' and ', ann*dble(pbc_111)/(dble(pbc_111)+1.d0)
        write(*,*) 'but no bigger than ...', x_high_limit/3.d0
        write(*,*) 'but no smaller than ...', x_low_limit/3.d0
        stop
      end if  
    end if  


 
    rcx(-2)=-dabs(dx2)
    rcy(-2)=-dabs(dx2)
    rcz(-2)=-dabs(dx2)

    rcx(-1)=-dabs(dx1)
    rcy(-1)=-dabs(dx1)
    rcz(-1)=-dabs(dx1)

    rcx(0)=0.d0
    rcy(0)=0.d0
    rcz(0)=0.d0

    rcx(1)=-rcx(-1)
    rcy(1)=-rcy(-1)
    rcz(1)=-rcz(-1)

    rcx(2)=-rcx(-2)
    rcy(2)=-rcy(-2)
    rcz(2)=-rcz(-2)
 
 return
end subroutine fill_the_111


 subroutine build_the_nsite_the_111(nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
 implicit none
  integer, intent(in) :: nsia,n1MAX,n2MAX,n3MAX
  integer, intent (in) :: ns1(nsia),ns2(nsia),ns3(nsia)
  integer  :: is,jj,n1pbc,n2pbc,n3pbc

  n1pbc=n1MAX+1
  n2pbc=n2MAX+1
  n3pbc=n3MAX+1


 do is=1,nsia
  do jj=-2,2,1
  ns1_111(jj,is)=ns1(is)+jj
  ns2_111(jj,is)=ns2(is)+jj
  ns3_111(jj,is)=ns3(is)+jj
  end do
 end do



 do is=1,nsia
  do jj=-2,2,1
   if (ns1_111(jj,is)< 0) ns1_111(jj,is)= ns1_111(jj,is)+n1pbc
   if (ns2_111(jj,is)< 0) ns2_111(jj,is)= ns2_111(jj,is)+n2pbc
   if (ns3_111(jj,is)< 0) ns3_111(jj,is)= ns3_111(jj,is)+n3pbc

   if (ns1_111(jj,is)>= n1pbc ) ns1_111(jj,is)= ns1_111(jj,is)-n1pbc
   if (ns2_111(jj,is)>= n2pbc ) ns2_111(jj,is)= ns2_111(jj,is)-n2pbc
   if (ns3_111(jj,is)>= n3pbc ) ns3_111(jj,is)= ns3_111(jj,is)-n3pbc
 end do
end do

 return

end subroutine  build_the_nsite_the_111



 subroutine build_the_nsite_the_111_migration(nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
 use migration_sia, ONLY : type_migration
  implicit none
  integer, intent(in) :: nsia,n1MAX,n2MAX,n3MAX
  integer, intent (in) :: ns1(nsia),ns2(nsia),ns3(nsia)
  integer  :: is,jj,n1pbc,n2pbc,n3pbc




  n1pbc=n1MAX+1
  n2pbc=n2MAX+1
  n3pbc=n3MAX+1


 do is=1,nsia
 if (type_migration(is)==0) then
  do jj=-2,2,1
   deb_ns1_111(jj,is)=ns1(is)+jj
   deb_ns2_111(jj,is)=ns2(is)+jj
   deb_ns3_111(jj,is)=ns3(is)+jj

   fin_ns1_111(jj,is)=ns1(is)+jj
   fin_ns2_111(jj,is)=ns2(is)+jj
   fin_ns3_111(jj,is)=ns3(is)+jj
  end do
 end if 

 if (type_migration(is)==1) then
  do jj=-2,2,1
   deb_ns1_111(jj,is)=ns1(is)+jj
   deb_ns2_111(jj,is)=ns2(is)+jj
   deb_ns3_111(jj,is)=ns3(is)+jj

   fin_ns1_111(jj,is)=ns1(is)+jj+1
   fin_ns2_111(jj,is)=ns2(is)+jj+1
   fin_ns3_111(jj,is)=ns3(is)+jj+1
   
  end do
 end if 

 end do

 write(*,*) 


 do is=1,nsia
  do jj=-2,2,1
   if (deb_ns1_111(jj,is)< 0) deb_ns1_111(jj,is)= deb_ns1_111(jj,is)+n1pbc
   if (deb_ns2_111(jj,is)< 0) deb_ns2_111(jj,is)= deb_ns2_111(jj,is)+n2pbc
   if (deb_ns3_111(jj,is)< 0) deb_ns3_111(jj,is)= deb_ns3_111(jj,is)+n3pbc

   if (deb_ns1_111(jj,is)>= n1pbc ) deb_ns1_111(jj,is)= deb_ns1_111(jj,is)-n1pbc
   if (deb_ns2_111(jj,is)>= n2pbc ) deb_ns2_111(jj,is)= deb_ns2_111(jj,is)-n2pbc
   if (deb_ns3_111(jj,is)>= n3pbc ) deb_ns3_111(jj,is)= deb_ns3_111(jj,is)-n3pbc

   if (fin_ns1_111(jj,is)< 0) fin_ns1_111(jj,is)= fin_ns1_111(jj,is)+n1pbc
   if (fin_ns2_111(jj,is)< 0) fin_ns2_111(jj,is)= fin_ns2_111(jj,is)+n2pbc
   if (fin_ns3_111(jj,is)< 0) fin_ns3_111(jj,is)= fin_ns3_111(jj,is)+n3pbc

   if (fin_ns1_111(jj,is)>= n1pbc ) fin_ns1_111(jj,is)= fin_ns1_111(jj,is)-n1pbc
   if (fin_ns2_111(jj,is)>= n2pbc ) fin_ns2_111(jj,is)= fin_ns2_111(jj,is)-n2pbc
   if (fin_ns3_111(jj,is)>= n3pbc ) fin_ns3_111(jj,is)= fin_ns3_111(jj,is)-n3pbc
 end do
end do

 return

end subroutine  build_the_nsite_the_111_migration




subroutine distord_lattice_111 (ntemp,nsia,rx,ry,rz)
implicit none
integer , intent(in) :: nsia,ntemp
real(8), intent(inout) :: rx(ntemp),ry(ntemp),rz(ntemp)
integer :: jj,is,i

 do jj=-2,2,1
  do is=1,nsia
     i=nsite_111(jj,is)
     rx(i) = rx(i)+rcx(jj)
     ry(i) = ry(i)+rcy(jj)
     rz(i) = rz(i)+rcz(jj)
   !  write(*,*) i, rcx(jj),rcy(jj),rcz(jj)
   end do
  end do

 return
end subroutine distord_lattice_111 

subroutine deb_distord_lattice_111_migration (ntemp,nsia,rxb,ryb,rzb,rx,ry,rz)
implicit none
integer , intent(in) :: nsia,ntemp
real(8), intent(in) :: rxb(ntemp),ryb(ntemp),rzb(ntemp)
real(8), intent(out):: rx(ntemp),ry(ntemp),rz(ntemp)
 
integer :: jj,is,i

 do jj=-2,2,1
  do is=1,nsia
     i=deb_nsite_111(jj,is)
     rx(i) = rxb(i)+rcx(jj)
     ry(i) = ryb(i)+rcy(jj)
     rz(i) = rzb(i)+rcz(jj)
  !   write(*,*) i, rcx(jj),rcy(jj),rcz(jj)
   end do
  end do


            write(*,'("deb side ",3f15.5)') rxb(7),ryb(7),rzb(7)


 return
end subroutine deb_distord_lattice_111_migration 

subroutine fin_distord_lattice_111_migration (ntemp,nsia,rxb,ryb,rzb,rx,ry,rz)
implicit none
integer , intent(in) :: nsia,ntemp
real(8), intent(in) :: rxb(ntemp),ryb(ntemp),rzb(ntemp)
real(8), intent(out) :: rx(ntemp),ry(ntemp),rz(ntemp)
 
integer :: jj,is,i

 do jj=-2,2,1
  do is=1,nsia
     i=fin_nsite_111(jj,is)
     rx(i) = rxb(i)+rcx(jj)
     ry(i) = ryb(i)+rcy(jj)
     rz(i) = rzb(i)+rcz(jj)


  !   write(*,*) i, rcx(jj),rcy(jj),rcz(jj)
   end do
  end do


 return
end subroutine fin_distord_lattice_111_migration 
end module the_111_smooth
