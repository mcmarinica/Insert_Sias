!=======================================================================
! kazuto_style.f90
!
! Type 67: a box with edges along [110], [-110] and [001] cut out of a
! crystal that holds a C15 cluster.
!
!   define_center                     corners of the cut in the (001) plane
!   surface_of_three_points_in_plane  area of a triangle
!   solid_and_defect_cut              marks the atoms inside the cut
!   write_out_kazuto                  builds the cut and writes the k_*.gin files
!=======================================================================

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

!=======================================================================

subroutine  surface_of_three_points_in_plane(xA,yA,xB,yB,xC,yC,surface_ij)
!object:
!       return the surface of a triangle having the theree points A,B,C in the corners.
!input:
!       the 2D coordinates of the point A: xA  yA
!       the 2D coordinates of the point B: xB  yB
!       the 2D coordinates of the point C: xC  yC
!output:
!       the surface of triangle ABC



implicit none
real(8),intent(in)   :: xA,yA,xB,yB,xC,yC
real(8), intent(out) :: surface_ij

!and the survivors _formula:

surface_ij=dabs(0.5d0*((xA-xC)*(yB-yA)-(xA-xB)*(yC-yA)))

return

end subroutine surface_of_three_points_in_plane

!=======================================================================

subroutine solid_and_defect_cut (a0, surface, new_solid, gndatoms, nsize, &
             rx,ry,rz, isite_surface, is_tot, izcoordinate, plane_corner)
!----------------------------------------------------------------------------
!object:  mark the atoms which are inside in the Kazuto's cut paralelipiped
!method:
!         the point P is  projected in (XY) plane
!         the sum of surfaces of triangles SurfaceP=Sum_i {Corner_i, Corner_i+1, Pprojection} is computed
!         if the:
!                   (Surface_base = SurfaceP) is true
!                   (Zcoordinate_of_P is between the Z of paralelipiped)  is true
!         then the point P is inside the paralelipiped
!input:
!output:
!          isite_surface(i)=0 not in paralelipiped
!          isite_surface(i)=i is included in paralelipiped
!-----------------------------------------------------------------------------
implicit none

integer, intent(in)  :: new_solid, gndatoms,nsize,izcoordinate
real(8), intent(in)  :: a0, surface
real(8), dimension(3,4), intent(in) :: plane_corner
real(8), dimension(nsize), intent(in) :: rx,ry,rz
integer, dimension(nsize),intent(out) :: isite_surface
integer, intent(out)  :: is_tot
real(8) :: surface_test, surface_kl
integer ::  itest_c15,k ,l,i


            is_tot=0
            itest_c15=0
            isite_surface(:)=0

             do i=1,new_solid
                surface_test=0.d0
                 do k=1,4
                  l=k+1
                   if (l==5) l=1
                      call surface_of_three_points_in_plane(       rx(i),ry(i),     &
                                    a0*dble(plane_corner(1,k))/2.d0,a0*dble(plane_corner(2,k))/2.d0,      &
                                    a0*dble(plane_corner(1,l))/2.d0,a0*dble(plane_corner(2,l))/2.d0,surface_kl)
                      surface_test=surface_test+surface_kl

                  !write(*,*) surface_test/(a0/2.d0)**2, surface_kl/(a0/2.d0)**2
                 end do

                 if (dabs(surface_test/(a0/2.d0)**2-surface)<=0.1d0) then
                  if ((dble(izcoordinate)+0.1d0)>=rz(i)/(a0/2.d0)) then
                    isite_surface(i)=1
                    is_tot=is_tot+1
                    if (i<=gndatoms) then
                      itest_c15=itest_c15+1
                    end if
                !debug   write(26,*)  rx(i), ry(i)
                  end if
                 end if
             !debug write (27,*) rx(i), ry(i)

             !write(*,*)i, surface, surface_test/(a0/2.d0)**2
             end do    ! i=new_solid
             if (itest_C15/=gndatoms) then
               write(*,*) 'WARNING NOT ALL THE C15 atoms are INCLUDED IN THE CUT'
             end if


return
end subroutine solid_and_defect_cut

!=======================================================================

subroutine write_out_kazuto (a0, n1MAX, n3MAX, ntemp, new_solid, gndatoms, nsize, nsia, nsite, &
                             rx, ry, rz, rxb, ryb, rzb)
!----------------------------------------------------------------------------
!object:  cut a box with edges along [110], [-110] and [001] out of the crystal
!         that holds the C15 cluster, and the same box out of the reference crystal
!input:
!          a0                lattice parameter
!          n1MAX, n3MAX      largest integer coordinates of the crystal along x and z
!          ntemp             atoms of the reference crystal
!          new_solid         atoms of the crystal with the cluster
!          gndatoms          inserted atoms; they come first in rx,ry,rz
!          nsite             atom at the centre of each C15 unit
!          rx,ry,rz          crystal with the cluster
!          rxb,ryb,rzb       reference crystal
!output:
!          k_dfct.gin        the cut of the crystal with the cluster
!          k_bulk.gin        the cut of the reference crystal
!          k_dfct_vide.gin   the inserted atoms and the centres of the units
!
!          The size of the cut is the constant pdist below, in units of a0/2.
!-----------------------------------------------------------------------------
use math
implicit none

integer, intent(in)  :: n1MAX, n3MAX, ntemp, new_solid, gndatoms, nsize, nsia
integer, intent(in)  :: nsite(nsia)
real(8), intent(in)  :: a0
real(8), dimension(nsize), intent(in) :: rx, ry, rz, rxb, ryb, rzb

integer :: pdist, izcoordinate, is_tot, is_totb, i, j, is, it_a, i_count
integer, dimension(:), allocatable :: isite_surface, isite_surfaceb
real(8), dimension(3,4) :: plane_corner
real(8) :: surface
real(8), dimension(3)   :: shift, temp, ctemp
real(8), dimension(3,3) :: amat, inv_amat, abase, inv_abase
real(8) :: dist_current, dist_ini

  allocate (isite_surface(nsize), isite_surfaceb(nsize))

  !pdist should be 2*k number is we want a complete cubic cell
  !pdist=(n1MAX-1)/2
  !pdist=48  !10 nm with 48-ABAB 49-ABA ! with cube.x = 89a0
  !pdist=18  !3  nm ! with cube.x = 51a0 or 49a0 or 53a0
  pdist=25    !5 nm !    24-ABAB  25-ABA   with cube.x = 89a0
  !pdist=5
  !if pdist is 2*k+1 the cell is incomplete ...

   !
   if ((pdist>(n1MAX-1)/2).or.(pdist<=0))  then
     write (*,'(" p should be < (n1MAX-1)/2; and > 0...:",2i7)') pdist,(n1MAX-1)/2
     write (*,*) 'stop in put_vct'
     stop
   end if
   !
   if (mod(pdist,2)/=0) then
     write (*,'(" p should be 2*k, k is an integer.....: ",2i7)') pdist,(n1MAX-1)/2
     !write (*,*) 'stop in put_vct'
     !stop
   end if

   call define_center (pdist,n1MAX,surface,plane_corner)
   izcoordinate=nint(dble(n1MAX-1)/2.d0*dsqrt(2.d0))
   if (mod(izcoordinate,2)/=0) then
     izcoordinate=izcoordinate+1
   end if
   izcoordinate = izcoordinate+1
   if (izcoordinate.gt.(n3MAX-1)) then
    write(*,*) 'At least the n3MAX should be increased ....', n3MAX
     stop
   end if

  call solid_and_defect_cut (a0, surface, new_solid, gndatoms, nsize, rx,ry,rz, isite_surface, &
                             is_tot,izcoordinate, plane_corner)
  call solid_and_defect_cut (a0, surface, ntemp, gndatoms, nsize, rxb,ryb,rzb, isite_surfaceb, &
                             is_totb,izcoordinate, plane_corner)
   open (unit=44, file='k_bulk.gin',status='unknown')
   open (unit=45, file='k_dfct.gin',status='unknown')
   open (unit=46, file='k_dfct_vide.gin',status='unknown')
   write(46,*) '1 1 1'
   write(45,*) '1 1 1'
   write(44,*) '1 1 1'

   shift(:)=(/ 1.d0,1.d0,0.d0 /)
   abase(:,1) = (a0/2.d0)*(plane_corner(:,2)-plane_corner(:,1) + shift(:))
   shift(:) = (/ -1.d0,1.d0,0.d0 /)
   abase(:,2) = (a0/2.d0)*(plane_corner(:,4)-plane_corner(:,1) + shift(:))
   abase(:,3) =  (a0/2.d0)*(/ 0.d0, 0.d0, dble(izcoordinate)+1.d0 /)

   abase(:,1) = (/ dsqrt(DOT_PRODUCT(abase(:,1),abase(:,1))), 0.d0, 0.d0 /)
   abase(:,2) = (/ 0.d0, dsqrt(DOT_PRODUCT(abase(:,2),abase(:,2))),       0.d0 /)
   abase(:,3) = (/ 0.d0, 0.d0, dsqrt(DOT_PRODUCT(abase(:,3),abase(:,3)))  /)

   call MatInv(abase,inv_abase)
   write(45,'(3f15.8)') abase(:,1)
   write(45,'(3f15.8)') abase(:,2)
   write(45,'(3f15.8)') abase(:,3)
   write(45,'(i6)') is_tot

   write(44,'(3f15.8)') abase(:,1)
   write(44,'(3f15.8)') abase(:,2)
   write(44,'(3f15.8)') abase(:,3)
   write(44,'(i6)') is_totb

   write(46,'(3f15.8)') abase(:,1)
   write(46,'(3f15.8)') abase(:,2)
   write(46,'(3f15.8)') abase(:,3)
   write(46,'(i6)') gndatoms+nsia

   amat(1,:) = (/  1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0),  0.d0 /)
   amat(2,:) = (/ -1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0),  0.d0 /)
   amat(3,:) = (/             0.d0,            0.d0,  1.d0 /)

   call MatInv(amat,inv_amat)

   shift(:) = plane_corner(:,1)*(a0/2.d0)

   !The suggested C15 center:
    do j=1,3
     temp(j)=shift(j) + amat(1,j)*abase(1,1)/2.d0+ amat(2,j)*abase(2,2)/2.d0+amat(3,j)*abase(3,3)/2.d0
    end do
    write(*,'("The suggested center for C15.....:",3i6)') int(temp(:)/(a0/2.d0))

   do i=1,gndatoms
        temp(:) = (/ rx(i), ry(i) , rz(i) /)
        !transfer into the new  oXYZ
        do j=1,3
         ctemp(j) =DOT_PRODUCT(inv_amat(:,j),temp(:)-shift(:))
        end do

        if (isite_surface(i)==1) then
        !    write(55,'(i6,a3,3f10.4,2f7.2)')i,' Fe',ctemp(1:3),1.d0,0.5d0
          write(45,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
          write(46,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
        end if
  end do

  do is=1,nsia
        i=nsite(is)
        temp(:) = (/ rxb(i), ryb(i) , rzb(i) /)
        !transfer into the new  oXYZ
        do j=1,3
         ctemp(j) =DOT_PRODUCT(inv_amat(:,j),temp(:)-shift(:))
        end do

        !    write(55,'(i6,a3,3f10.4,2f7.2)')i,' Fe',ctemp(1:3),1.d0,0.5d0
          write(45,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
          write(46,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
  end do

   i_count=0
   do i=gndatoms+1,new_solid

     dist_ini=99999999.9
      !
      do is=1,nsia
        it_a=nsite(is)
        dist_current=dsqrt((rx(i)-rxb(it_a))**2+(ry(i)-ryb(it_a))**2+(rz(i)-rzb(it_a))**2)
        if (dist_current < dist_ini) dist_ini=dist_current
      end do
      !

     if (dist_ini  > 0.01) then
        i_count=i_count+1
        temp(:) = (/ rx(i), ry(i) , rz(i) /)
        !transfer into the new  oXYZ
        do j=1,3
         ctemp(j) =DOT_PRODUCT(inv_amat(:,j),temp(:)-shift(:))
        end do

        if (isite_surface(i)==1) then
        !    write(55,'(i6,a3,3f10.4,2f7.2)')i,' Fe',ctemp(1:3),1.d0,0.5d0
          write(45,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
        end if
      end if

  end do

   do i=1,ntemp
      temp(:) = (/ rxb(i), ryb(i) , rzb(i) /)
      !transfer into the new reper
      do j=1,3
        ctemp(j) =DOT_PRODUCT(inv_amat(:,j),temp(:)-shift(:))
      end do
      !
      if (isite_surfaceb(i)==1) then
          write(44,'(a3,3(f18.13),"   1")') '   ',  MatMul(inv_abase(:,:), ctemp(:))
      end if
  end do
  close (44)
  close (45)

return
end subroutine write_out_kazuto

!=======================================================================
