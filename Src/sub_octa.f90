!#####################################################################


!#####################################################################
   subroutine fix_ico_type (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
   use icosaedre, only:  type_icos
   integer  :: nsia, n1MAX,n2MAX,n3MAX,i
   integer  :: ns1(nsia),ns2(nsia),ns3(nsia), i1,i2,i3

   if (allocated(type_icos)) deallocate(type_icos);  allocate (type_icos(nsia))
      do i=1,nsia

        i1 =  mod(ns1(i),2)
        i2 =  mod(ns2(i),2)
        i3 =  mod(ns3(i),2)

        if (i1==0) type_icos(i)=2
        if (i1==1) type_icos(i)=1

        if ((i1 /= i2) .or. (i1/= i3)) then
            write(*,*) 'sia site ', i, 'all the line should have the same parity'
            write(*,*) 'ns1, ns2, ns3', ns1(i),ns2(i), ns3(i)
            stop
        end if
      end do
   return
   end


   subroutine fix_ico_type_mixed (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
   use icosaedre, only:  type_icos, nsia_ico, nsia_dumbell
   integer  :: nsia, n1MAX,n2MAX,n3MAX,i
   integer  :: ns1(nsia),ns2(nsia),ns3(nsia), i1,i2,i3

   if (allocated(type_icos)) deallocate(type_icos);  allocate (type_icos(nsia))
      do i=1,nsia_ico

        i1 =  mod(ns1(i),2)
        i2 =  mod(ns2(i),2)
        i3 =  mod(ns3(i),2)

        if (i1==0) type_icos(i)=2
        if (i1==1) type_icos(i)=1

        if ((i1 /= i2) .or. (i1/= i3)) then
            write(*,*) 'sia site ', i, 'all the line should have the same parity'
            write(*,*) 'ns1, ns2, ns3', ns1(i),ns2(i), ns3(i)
            stop
        end if
      end do


   return
   end




  subroutine read_config_sia_octa (nsia,type_octa,  &
                    ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)

  integer  :: nsia, n1MAX,n2MAX,n3MAX,i
  integer  :: ns1(nsia),ns2(nsia),ns3(nsia),type_octa(nsia)


!          do i=1,nsia
!           read (11,*) ns1(i),ns2(i),ns3(i)
!           write(*,*) 'ns22',ns1(i),ns2(i),ns3(i)
!            if ( (ns1(i) > n1MAX) .or.(ns2(i) > n2MAX) .or. &
!              (ns3(i) > n3MAX) ) then
!            print*, 'Big problem in read_config_sia'
!            print*, 'The cell is not big enough. Increase the size!'
!              write(*,*) '<read_config_sia_octa>'
!            stop
!            end if
!          end do
     do i=1,nsia
      read(11,*) type_octa(i)
      !debug write(*,*) 'ns22 HERRRRRRRRRRRE',  type_octa(i)
      if ( (type_octa(i).lt.1).or.((type_octa(i).gt.4)) ) then
      write(*,*) 'There is no implementation for this octa type', &
        type_octa(i)
       write(*,*) '<read_config_sia_octa>'
       stop
      end if
     end do

  return
  end


  subroutine read_config_sia_ico (nsia,nd1,nd2,nd3)
  use icosaedre, only : type_dumb
  integer  :: nsia, i
  integer  :: nd1(nsia),nd2(nsia),nd3(nsia)


     do i=1,nsia
      read(11,*) nd1(i),nd2(i), nd3(i), type_dumb(i)
      !debug write(*,*) 'ns22 HERRRRRRRRRRRE',  type_octa(i)
      if ( (type_dumb(i).lt.1).or.((type_dumb(i).gt.4)) ) then
      write(*,*) 'There is no implementation for this octa type'
      write(*,*) '<read_config_sia_ico>'
       stop
      end if
     end do

  return
  end



!#####################################################################

  subroutine  lattice_compute_different_removeatoms (nsia, &
      nlac_remove, nsite_octa_vac_remove,indx,list_no_order,ndatoms)


  implicit none
    !intent IN/OUT
  integer, intent(in) :: nsia,nlac_remove
    integer, intent(out):: ndatoms
    integer, intent(out) :: indx(nlac_remove)
  real(8), intent(out)  :: list_no_order(nlac_remove)
  !local variables
    integer :: ii,i_double,icnt,i
  integer :: nsite_octa_vac_remove(10*nsia)

  if (nlac_remove.gt.10*nsia) then
   write(*,*) 'nlac_remove is greater than 10*nsia'
   write(*,*) 'Problems in dimension of the vector nsite_octa_vac_remove !!'
   stop
  end if

    icnt=0
    do ii=1,nlac_remove
      list_no_order(ii)=nsite_octa_vac_remove(ii)
      !debug write(*,*) 'list', list_no_order(ii),ii, nlac_remove, nlac_remove
    end do

   call indexx (nlac_remove,list_no_order,indx)


   i_double=0
   !debug write(*,*) 'list', list_no_order(indx(1)),i_double
   do    i=1,nlac_remove-1
     if (dabs(list_no_order(indx(i+1))-list_no_order(indx(i))).lt.0.01d0) then
       i_double = i_double+1
     end if
     !debug write(*,*) 'list', list_no_order(indx(i+1)),i_double
   end do

    do ii=1,nlac_remove
       if(dabs(list_no_order(ii)-nsite_octa_vac_remove(ii)).gt.1d-5) then
         write(*,*) '!!!! in list_no_order and nsite_octa_vac_remove'
       end if
    end do


     write(*,*) 'The number of total vacancies (which can be superposed).....: ',nlac_remove
     write(*,*) 'The number of DIFFERENT VACANCIES left in the input ........: ',nlac_remove-i_double
     ndatoms=nlac_remove-i_double

   return
   end

!#####################################################################

 subroutine  lattice_compute_atoms_to_remove  (atoms_to_remove,&
               atoms_to_remove_type,ndatoms,nlac_remove,         &
                list_no_order,nsite_octa_vac_remove_type,indx)

    integer, intent(in)  :: ndatoms,nlac_remove
  real(8), intent(in)  :: list_no_order(nlac_remove)
  integer, intent(in)  :: indx(nlac_remove),nsite_octa_vac_remove_type(nlac_remove)
  integer, intent(out)  :: atoms_to_remove_type(ndatoms)
  real(8), intent(out)   :: atoms_to_remove(ndatoms)
    integer :: i_double,icnt,i


  atoms_to_remove(1)=list_no_order(indx(1))
  atoms_to_remove_type(1)=nsite_octa_vac_remove_type(indx(1))
  i_double=0
  icnt=1
  do    i=1,nlac_remove-1
    if (dabs(list_no_order(indx(i+1))-list_no_order(indx(i))).lt.0.01d0) then
      i_double = i_double+1
    else
      icnt=icnt+1
      atoms_to_remove(icnt)=list_no_order(indx(i+1))
      atoms_to_remove_type(icnt)=nsite_octa_vac_remove_type(indx(i+1))
    end if
  end do

  !debug do i=1,icnt
  !debug    write(*,*) 'toremove', i, atoms_to_remove(i), ndatoms
  !debug end do

  !debug write(*,*) 'list', icnt,i_double,10*nsia-i_double


  return
  end
!#####################################################################

  subroutine  offlattice_get_coordinates_and_order(a0,rxgaos,rygaos,rzgaos,&
                                                   rx,    ry,    rz,    &
                 offlattice_no_order,indx,nsite,type_octa,nsia,nsize,gndatoms,natoms_extra)

   use octahedre
   implicit none

   integer, intent(in) :: nsia,nsize
   integer, intent(in) :: nsite(nsia),type_octa(nsia)
   real(8), intent(in)   :: a0
   real(8), dimension(nsize)  :: rx,ry,rz
   integer,  intent(out) :: indx (12*nsia),gndatoms,natoms_extra

   real(8), dimension(12*nsia),  intent(out) :: offlattice_no_order,rxgaos,rygaos,rzgaos
   integer :: i,it_octo,i_double,icnt,ii,jj
   real(8)   :: a02,alpha,origin(3)


   a02=a0/2.d0
   alpha=0.5d0
   icnt=0
  do i=1,nsia
   it_octo=type_octa(i)
     origin(1:3)=(/ rx(nsite(i)),  ry(nsite(i)), rz(nsite(i)) /)
       if ( (it_octo.eq.1).or.(it_octo.eq.2) )then
      do ii=1,4  !loop over the 4 corners of the tetrahedron
       do jj=1,3 !loop over the three atoms located in the corner of the tetrahedron
           icnt=icnt+1
         rxgaos(icnt)=a02*(alpha*octo(it_octo)%gaos(1,jj,ii)+octo(it_octo)%atom(1,ii))+ origin(1)
         rygaos(icnt)=a02*(alpha*octo(it_octo)%gaos(2,jj,ii)+octo(it_octo)%atom(2,ii))+ origin(2)
         rzgaos(icnt)=a02*(alpha*octo(it_octo)%gaos(3,jj,ii)+octo(it_octo)%atom(3,ii))+ origin(3)
         offlattice_no_order(icnt)=rxgaos(icnt)**3+10*rygaos(icnt)**3+100*rzgaos(icnt)**3
       end do
      end do
      else if (it_octo.eq.4) then

      do ii=1,1  !loop over the 4 corners of the tetrahedron
       do jj=1,3 !loop over the three atoms located in the corner of the tetrahedron
           icnt=icnt+1
         rxgaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(1,jj,ii)+octo(it_octo-2)%atom(1,ii))+ origin(1)
         rygaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(2,jj,ii)+octo(it_octo-2)%atom(2,ii))+ origin(2)
         rzgaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(3,jj,ii)+octo(it_octo-2)%atom(3,ii))+ origin(3)
         offlattice_no_order(icnt)=rxgaos(icnt)**3+10*rygaos(icnt)**3+100*rzgaos(icnt)**3
       end do
      end do

      do ii=3,3  !loop over the 4 corners of the tetrahedron
       do jj=2,2 !loop over the three atoms located in the corner of the tetrahedron
           icnt=icnt+1
         rxgaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(1,jj,ii)+octo(it_octo-2)%atom(1,ii))+ origin(1)
         rygaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(2,jj,ii)+octo(it_octo-2)%atom(2,ii))+ origin(2)
         rzgaos(icnt)=a02*(alpha*octo(it_octo-2)%gaos(3,jj,ii)+octo(it_octo-2)%atom(3,ii))+ origin(3)
         offlattice_no_order(icnt)=rxgaos(icnt)**3+10*rygaos(icnt)**3+100*rzgaos(icnt)**3
       end do
      end do

      end if
  end do

   natoms_extra=icnt
!old       call indexx (12*nsia, offlattice_no_order,indx)
   call indexx (icnt, offlattice_no_order,indx)

   i_double=0
   !debug write(*,'("listgaos",f14.3,i4,3f12.3)') &
   !debug offlattice_no_order(indx(1)),i_double,rxgaos(1),rygaos(1),rzgaos(1)
   do    i=1,icnt-1
     !old was 1.1
   if (dabs(offlattice_no_order(indx(i+1))-offlattice_no_order(indx(i))).lt.1.1d0) then
   i_double = i_double+1
   end if
   !debug write(*,'("listgaos",f14.3,i4,3f12.3)') &
   !debug offlattice_no_order(indx(i+1)),i_double,&
   !debug rxgaos(indx(i+1)),rygaos(indx(i+1)),rzgaos(indx(i+1))
   end do


     write(*,*) 'The number of DIFFERENT GAOS atoms which must be removed ....: ',i_double
     write(*,*) 'The number of DIFFERENT GAOS atoms left in input ............: ',icnt-i_double
     gndatoms=icnt-i_double

   return
   end


!#####################################################################

subroutine  ico_offlattice_get_coordinates_and_order(a0,rxgaos,rygaos,rzgaos,&
                                                   rx,    ry,    rz,ns1,ns2,ns3,    &
                 offlattice_no_order,indx,nsite,type_octa,nsia,nsize,gndatoms,natoms_extra)

use icosaedre
implicit none

integer, intent(in) :: nsia,nsize
integer, intent(in) :: nsite(nsia),type_octa(nsia), ns1(nsia), ns2(nsia), ns3(nsia)
real(8), intent(in)   :: a0
real(8), dimension(nsize)  :: rx,ry,rz
integer,  intent(out) :: indx (13*nsia),gndatoms,natoms_extra

real(8), dimension(13*nsia),  intent(out) :: offlattice_no_order,rxgaos,rygaos,rzgaos
integer :: i,it_icos,i_double,icnt,ii,jj
real(8)   :: a02,alpha,origin(3)


a02=a0/2.d0
alpha=1.0d0
icnt=0
do i=1,nsia
it_icos=type_icos(i)
  origin(1:3)=(/ ns1(i),  ns2(i), ns3(i) /)
    do jj=1,12 !loop over the three atoms located in the corner of the tetrahedron
       icnt=icnt+1
       rxgaos(icnt)=a02*(alpha*icos(it_icos)%sia(1,jj)+ origin(1))
       rygaos(icnt)=a02*(alpha*icos(it_icos)%sia(2,jj)+ origin(2))
       rzgaos(icnt)=a02*(alpha*icos(it_icos)%sia(3,jj)+ origin(3))
       offlattice_no_order(icnt)=rxgaos(icnt)**3+10*rygaos(icnt)**3+100*rzgaos(icnt)**3
       write(*,'("ig", i6,f20.7,3f18.6)') icnt, offlattice_no_order(icnt), rxgaos(icnt), rygaos(icnt), rzgaos(icnt)
   end do
   !if (it_icos ==  1) then
     do jj=1,1 !loop over the three atoms located in the corner of the tetrahedron
       icnt=icnt+1
       rxgaos(icnt)=a02*(alpha*icos(it_icos)%add(1,jj)+ origin(1))
       rygaos(icnt)=a02*(alpha*icos(it_icos)%add(2,jj)+ origin(2))
       rzgaos(icnt)=a02*(alpha*icos(it_icos)%add(3,jj)+ origin(3))
       offlattice_no_order(icnt)=rxgaos(icnt)**3+10*rygaos(icnt)**3+100*rzgaos(icnt)**3
       write(*,'("ig", i6,f20.7,3f18.6)') icnt, offlattice_no_order(icnt), rxgaos(icnt), rygaos(icnt), rzgaos(icnt)
     end do
   !end if

end do

 natoms_extra=icnt
 if (natoms_extra.ne.13*nsia) then
   write(*,*) '<ico_offlattice_get_coordinates_and_order>  no of offlatice atoms is wrong', 13*nsia, natoms_extra
   stop
 end if

 call indexx (13*nsia, offlattice_no_order,indx)


 do jj=1,natoms_extra
   ii=indx(jj)
   write(*,'("gg", i6,f20.7,3f18.6)') ii, offlattice_no_order(ii), rxgaos(ii), rygaos(ii), rzgaos(ii)
 end do

 i_double=0
 !debug write(*,'("listgaos",f14.3,i4,3f12.3)') &
 !debug offlattice_no_order(indx(1)),i_double,rxgaos(1),rygaos(1),rzgaos(1)
 do    i=1,icnt-1
   !old was 1.1
 if (dabs(offlattice_no_order(indx(i+1))-offlattice_no_order(indx(i))).lt.1.1d0) then
 i_double = i_double+1
 end if
 !debug write(*,'("listgaos",f14.3,i4,3f12.3)') &
 !debug offlattice_no_order(indx(i+1)),i_double,&
 !debug rxgaos(indx(i+1)),rygaos(indx(i+1)),rzgaos(indx(i+1))
 end do


   write(*,*) 'The number of removed GAOS ..................................: ',i_double
   write(*,*) 'The number of DIFFERENT GAOS atoms left in input ............: ',icnt-i_double
   gndatoms=icnt-i_double

return
end subroutine  ico_offlattice_get_coordinates_and_order

subroutine ico_offlattice_compute_atoms_to_remove(gndatoms,rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos,&
                      offlattice_no_order,indx,nsia,natoms_extra)

implicit none

integer, intent(in) :: nsia
integer,  intent(in) :: indx (13*nsia),gndatoms,natoms_extra

real(8), dimension(13*nsia),    intent(in) :: offlattice_no_order,temprxgaos,temprygaos,temprzgaos
real(8), dimension(gndatoms),  intent(out) :: rxgaos,rygaos,rzgaos
integer :: i_double,i,icnt

i_double=0
icnt=1

 rxgaos(1)=temprxgaos(indx(1))
 rygaos(1)=temprygaos(indx(1))
 rzgaos(1)=temprzgaos(indx(1))
do    i=1,natoms_extra-1
  if (dabs(offlattice_no_order(indx(i+1))-offlattice_no_order(indx(i))).lt.0.1d0) then
     i_double = i_double+1
  else
    icnt=icnt+1
    rxgaos(icnt)=temprxgaos(indx(i+1))
    rygaos(icnt)=temprygaos(indx(i+1))
    rzgaos(icnt)=temprzgaos(indx(i+1))
  end if
end do

do i=1,icnt
write(*,*) 'gaos', i, rxgaos(i), rygaos(i),rzgaos(i)
end do
  write(*,'("II Final number of removed GAOS  .: ",i6)') i_double
return
end subroutine ico_offlattice_compute_atoms_to_remove



subroutine offlattice_compute_atoms_to_remove(gndatoms,rxgaos,rygaos,rzgaos,temprxgaos,temprygaos,temprzgaos,&
                      offlattice_no_order,indx,nsia,natoms_extra)

implicit none

integer, intent(in) :: nsia
integer,  intent(in) :: indx (12*nsia),gndatoms,natoms_extra

real(8), dimension(12*nsia),    intent(in) :: offlattice_no_order,temprxgaos,temprygaos,temprzgaos
real(8), dimension(gndatoms),  intent(out) :: rxgaos,rygaos,rzgaos
integer :: i_double,i,icnt

i_double=0
icnt=1

 rxgaos(1)=temprxgaos(indx(1))
 rygaos(1)=temprygaos(indx(1))
 rzgaos(1)=temprzgaos(indx(1))
do    i=1,natoms_extra-1
  if (dabs(offlattice_no_order(indx(i+1))-offlattice_no_order(indx(i))).lt.0.1d0) then
i_double = i_double+1
else
 icnt=icnt+1
 rxgaos(icnt)=temprxgaos(indx(i+1))
 rygaos(icnt)=temprygaos(indx(i+1))
 rzgaos(icnt)=temprzgaos(indx(i+1))

end if
end do

  write(*,'("II Final number of removed GAOS  .: ",i6)') i_double
return
end subroutine offlattice_compute_atoms_to_remove

!#####################################################################




subroutine write_coordinates (ntemp,nsize,new_solid,rxgaos,rygaos,rzgaos,gndatoms,atoms_to_remove,ndatoms,&
                        rx,ry,rz)

use octahedre
implicit none
integer,intent(in)  :: ntemp,nsize,ndatoms,gndatoms
integer, intent(out)  :: new_solid
double precision,dimension(gndatoms),intent(in)     :: rxgaos,rygaos,rzgaos
double precision,dimension(ndatoms),intent(in)     :: atoms_to_remove(ndatoms)
double precision,dimension(nsize),intent(inout)     :: rx,ry,rz
double precision,dimension(nsize)     :: rxtemp,rytemp,rztemp
integer :: i,icnt,is,icard,i_test_card
integer :: i1,i2
!
!
icnt=0
rxtemp(:)=99999999.d0
rytemp(:)=99999999.d0
rztemp(:)=99999999.d0
!ntemp   - the bulk reference
!ndatoms - the atoms that should be removed
!gndatoms- the extratoms as gaos, sias etc
!debug write(*,*) 'write_coordinates  ntemp, ndatoms, gndatoms :', ntemp, ndatoms, gndatoms

!writting the SIAs in 77 file
do i=1,gndatoms
  icnt=icnt+1
  rxtemp(icnt)=rxgaos(i)
  rytemp(icnt)=rygaos(i)
  rztemp(icnt)=rzgaos(i)
  write (77,'("atom Fe1   ",i6, 3f15.4)') icnt, rxgaos(icnt),rygaos(icnt),rzgaos(icnt)
end do


i1=0
i2=0
i_test_card=0
do i=1,ntemp
  icard=0
  do is=1,ndatoms
    if (i==nint(atoms_to_remove(is))) then
      i_test_card=i_test_card+1
      icard=1
      i1=i1+1
    endif
  end do
  if (icard==0) then
    i2=i2+1
    icnt=icnt+1
    rxtemp(icnt)=rx(i)
    rytemp(icnt)=ry(i)
    rztemp(icnt)=rz(i)
  end if
end do

  if (icnt.gt.nsize) then
   write(*,*) '<write_cordinates> Vector size not correct...icnt vs nsize:',icnt,nsize
   write(*,*) '<write_cordinates> icnt should be lower than nsize:',icnt,nsize
 ! stop
 end if

   write(*,*) i1,i2,gndatoms, i1+i2+gndatoms
   write(*,*) '<write_coordinates> Final number of atoms.:',icnt
 new_solid=icnt
 !in rx the first new_solid atoms are the real atoms of the final sistem. All others up to nsize are ghost atoms
 rx(1:nsize)=rxtemp(1:nsize)!*volc
 ry(1:nsize)=rytemp(1:nsize)!*volc
 rz(1:nsize)=rztemp(1:nsize)!*volc

 return
 end subroutine write_coordinates

subroutine dumbell_ico(rxb,ryb,rzb,nr, nd1,nd2,nd3,nd, a02)
use icosaedre, only : nsia_dumbell, r_extra_ico, nsite_dumb_remove, rdumb, type_dumb, nsia_dumbell_real_ico, nsia_center_real_ico
implicit none
integer :: nr, nd
real(8) :: rxb(nr), ryb(nr), rzb(nr),a02
integer :: nd1(nd), nd2(nd), nd3(nd)
integer :: i, ir, ii, icount

icount=0
do i=1,nsia_dumbell
   if (type_dumb(i)<=3) icount = icount + 2
   if (type_dumb(i)==4) icount = icount + 1
end do
   if (allocated(r_extra_ico)) deallocate(r_extra_ico) ; allocate(r_extra_ico(3,icount))
   if (allocated(nsite_dumb_remove)) deallocate(nsite_dumb_remove) ; allocate(nsite_dumb_remove(nsia_dumbell))



icount=0
do ir= 1, size(rxb)
      do ii=1,size(nd1)
        if (type_dumb(ii)/=4) then
        if ((nint(rxb(ir)/a02)==nd1(ii)).and.(nint(ryb(ir)/a02)==nd2(ii)).and.(nint(rzb(ir)/a02)==nd3(ii))) then
          nsite_dumb_remove(ii) = ir
          icount = icount +1
        end if
        end if
      end do
end do
nsia_dumbell_real_ico = icount
nsia_center_real_ico = size(nd1) - icount

write (*,*) nsia_dumbell_real_ico, nsia_center_real_ico

if (icount > size(nd1)) then
   write(*,*) 'HUge problem: icount > nd1. icount = the number of dumbbles , b   '
end if


icount=0
do i=1, nsia_dumbell
    if(type_dumb(i)<=3) then
      icount = icount+1
      r_extra_ico(1,icount) = rxb(nsite_dumb_remove(i)) + rdumb(1,1,type_dumb(i))
      r_extra_ico(2,icount) = ryb(nsite_dumb_remove(i)) + rdumb(1,2,type_dumb(i))
      r_extra_ico(3,icount) = rzb(nsite_dumb_remove(i)) + rdumb(1,3,type_dumb(i))


      icount = icount+1
      r_extra_ico(1,icount) = rxb(nsite_dumb_remove(i)) + rdumb(2,1,type_dumb(i))
      r_extra_ico(2,icount) = ryb(nsite_dumb_remove(i)) + rdumb(2,2,type_dumb(i))
      r_extra_ico(3,icount) = rzb(nsite_dumb_remove(i)) + rdumb(2,3,type_dumb(i))
 end if
end do

if (icount /= (2*nsia_dumbell_real_ico)) then
   write (*,*) 'Huge problem in dumbell_ico 1 ', icount, nsia_dumbell_real_ico
   stop
end if

do i=1,nsia_dumbell
    if (type_dumb(i)==4) then
      icount = icount+1
      r_extra_ico(1,icount) = nd1(i)*a02
      r_extra_ico(2,icount) = nd2(i)*a02
      r_extra_ico(3,icount) = nd3(i)*a02
    end if
end do


if ((icount - 2*nsia_dumbell_real_ico) /= nsia_center_real_ico) then
   write (*,*) 'Huge problem in dumbell_ico 2 '
   stop
end if

return
end subroutine dumbell_ico


!#####################################################################

       subroutine write_out_octa (nall,ntemp,new_solid,gndatoms,volc,&
                       a0, rx, ry, rz, rxb, ryb, rzb, a, ndatoms,    &
                       atoms_to_remove, atoms_to_remove_type,nsite,nsia)

      ! new_solid = the total atoms in the file
      ! ntemp  = the bulk atoms
      ! nall = real and imaginary atoms
      ! 1 to gndatoms  = the SIAs atoms in the case of C15
      ! gndatoms+1 to new_solid = all the others atoms ( the perfect bulk)
      ! ndatoms = the no of vacancies ( the atoms which are removed due to C15 phase)
       use math
       use pbc
       implicit none
       integer :: nsia
       integer :: nsite(nsia)
       integer  :: nall,ntemp,i,in,icnt,icnt2,gndatoms,new_solid,is
       double precision     :: rx(nall),ry(nall),rz(nall),a0,a02,CELLDM
       double precision, intent(in) :: rxb(nall),ryb(nall),rzb(nall)
       double precision     :: x,y,z, volc , delta
       double precision,dimension(3,3) :: a
       integer, intent(in)            :: ndatoms
       double precision,intent(in)     :: atoms_to_remove(ndatoms)
       integer,intent(in)     :: atoms_to_remove_type(ndatoms)
       double precision :: xmin,ymin,zmin,xmax,ymax,zmax
       logical   :: lremove
       real(kind(1.d0)), dimension(3) :: xp,xc
       real(kind(1.d0))   :: dist_current, dist_ini
       integer :: it_a,i_count
       integer :: f_jsia_xyz, f_sia_xyz, f_usia_xyz, f_jusia_xyz, f_uvsia_xyz, f_uvbsia_xyz, &
                  f_uvbcsia_xyz, f_uvbcsia_sim_xyz, f_unit_cell,  f_geom_new_data, f_geom_org_data, &
                  f_new_ndm_scl,f_formd_new_str
       write(*,*) ' nall,ntemp,new_solid,gndatoms, ndatoms',nall,ntemp,new_solid,gndatoms, ndatoms

       CELLDM=1.0
       in=1
       a02=a0/2.d0



       f_jsia_xyz=11
       open(f_jsia_xyz,file='jsia.xyz')
       f_sia_xyz=12
       open(12,file='sia.xyz')

       f_usia_xyz=19
       open(f_usia_xyz,file='usia.xyz')

       f_jusia_xyz=16
       open(f_jusia_xyz,file='jusia.xyz')

       f_uvsia_xyz=20
       open(f_uvsia_xyz,file='uvsia.xyz')


       f_uvbsia_xyz=21
       open(f_uvbsia_xyz,file='uvbsia.xyz')

       f_uvbcsia_xyz=23
       open(f_uvbcsia_xyz,file='uvbcsia.xyz')


       f_uvbcsia_sim_xyz=24
       open(f_uvbcsia_sim_xyz,file='uvbcsia_sim.xyz')

       f_unit_cell=13
       open(f_unit_cell,file='unit_cell',form='formatted')

       f_geom_new_data=17
       open (f_geom_new_data,file='geom_new.data',form='formatted',status='unknown')


       f_geom_org_data=22
       open (f_geom_org_data,file='geom_org.data',form='formatted',status='unknown')

       f_new_ndm_scl=18
       open (f_new_ndm_scl,file='formd_new_ndm_scl.str',form='formatted',&
          status='unknown')

      f_formd_new_str=15
      open (f_formd_new_str,file='formd_new.str',form='formatted',status='unknown')

      write(f_formd_new_str,'(i8)') new_solid
      write(f_formd_new_str,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
      write(f_formd_new_str,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
      write(f_formd_new_str,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


      write(f_geom_new_data,'(i8)') new_solid
      write(f_geom_new_data,*) 1.0
      write(f_geom_new_data,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
      write(f_geom_new_data,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
      write(f_geom_new_data,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

      write(f_geom_org_data,'(i8)') ntemp
      write(f_geom_org_data,*) 1.0
      write(f_geom_org_data,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
      write(f_geom_org_data,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
      write(f_geom_org_data,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc


   write(42,'("1   1   1")')
   write(42,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
   write(42,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
   write(42,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
   write(42,'(i8)') ntemp


      write(f_new_ndm_scl,'("1   1   1")')
      write(f_new_ndm_scl,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
      write(f_new_ndm_scl,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
      write(f_new_ndm_scl,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
      write(f_new_ndm_scl,'(i8)') new_solid



10       format(I8,1X,D24.17,1X,D24.17,1X,D24.17,1X,I1,1X,I1,1X,I1)
          write(f_unit_cell,*) 'number of atoms in the unit cell '
          write(f_unit_cell,*) new_solid
          write(f_unit_cell,*) 'CELLDM= '
          write(f_unit_cell,*)  CELLDM
        write(f_unit_cell,*) 'position of the atoms in the unit cell'



       write(f_sia_xyz,'(i8)') new_solid
       write(f_sia_xyz,'(a)') ' '
       write(f_jsia_xyz,'(i8)') new_solid
       write(f_jsia_xyz,'(a)') ' '
       write(f_usia_xyz,'(i8)') gndatoms
       write(f_usia_xyz,'(a)') ' '
       write(f_jusia_xyz,'(i8)') gndatoms
       write(f_jusia_xyz,'(a)') ' '
       write(f_uvsia_xyz,'(i8)') gndatoms+ndatoms
       write(f_uvsia_xyz,'(a)') ' '

       write(f_uvbsia_xyz,'(i8)') gndatoms+ndatoms+ntemp
       write(f_uvbsia_xyz,'(a)') ' '



       xmin=1.d09
       ymin=1.d09
       zmin=1.d09

       xmax=-1.d09
       ymax=-1.d09
       zmax=-1.d09

       do i=1,ntemp
         x = rxb(i)*volc
         y = ryb(i)*volc
         z = rzb(i)*volc
         do is=1,ndatoms
           if (i==nint(atoms_to_remove(is))) then
             if (x .ge. xmax) xmax=x
             if (y .ge. ymax) ymax=y
             if (z .ge. zmax) zmax=z
             if (x .le. xmin) xmin=x
             if (y .le. ymin) ymin=y
             if (z .le. zmin) zmin=z
           endif
         end do
       end do


       icnt2=0
       delta=+0.2
       xmin=xmin-delta
       ymin=ymin-delta
       zmin=zmin-delta
       xmax=xmax+delta
       ymax=ymax+delta
       zmax=zmax+delta

       do i=1,ntemp
         x = rxb(i)*volc
         y = ryb(i)*volc
         z = rzb(i)*volc
         if (x.gt.xmin) then
          if (y.gt.ymin) then
           if (z.gt.zmin) then

            if (x.lt.xmax) then
             if (y.lt.ymax) then
              if (z.lt.zmax) then
               icnt2=icnt2+1
              end if
             end if
            end if

           end if
          end if
         end if
       end do

       write(f_uvbcsia_xyz,'(i8)') gndatoms+ndatoms+icnt2
       write(f_uvbcsia_sim_xyz,'(i8)') gndatoms+ndatoms+icnt2
       write(f_uvbcsia_xyz,'(a)') ' '
       !gndatoms - atoms in interstitial position - label  as Fe

       icnt=0
       do i=1,gndatoms
       icnt =icnt+1
         x = rx(i)*volc
         y = ry(i)*volc
         z = rz(i)*volc
        write(f_sia_xyz,'("Fe  ",3f18.13)') x,y,z
        write(f_jsia_xyz,'("Fe  ",3f18.13)') x,y,z
        write(f_usia_xyz,'("Fe  ",3f18.13)') x,y,z
        write(f_jusia_xyz,'("Fe  ",3f18.13)') x,y,z
        write(f_uvsia_xyz,'("Fe  ",3f18.13)') x,y,z
        write(f_uvbsia_xyz,'("Fe  ",3f18.13)') x,y,z
        write(f_uvbcsia_xyz,'("Fe  ",3f18.13)') x,y,z
        write(f_uvbcsia_sim_xyz,'("Fe  ",3f18.13)') x,y,z
        write(f_unit_cell,10) icnt,x,y,z,in,in,in
        write(f_formd_new_str,'(a3,3(f18.13))') 'Fe ', x,y,z
        write(f_geom_new_data,'(a3,3(f18.13))') '   ', x,y,z
        xp(1:3)=(/  x, y, z/)
        xc(1:3) = MatMul(inv_at(1:3,1:3)/volc, xp(1:3))
        write(f_new_ndm_scl,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
       end do




       do i=1,ntemp
             x = rxb(i)*volc
             y = ryb(i)*volc
             z = rzb(i)*volc
            write(f_uvbsia_xyz,'("Cu  ",3f18.13)') x,y,z
        do is=1,ndatoms
         if (i==nint(atoms_to_remove(is))) then

           if (atoms_to_remove_type(is).eq.70) then
            write(f_uvsia_xyz,'("He  ",3f18.13)') x,y,z
            write(f_uvbsia_xyz,'("He  ",3f18.13)') x,y,z
            write(f_uvbcsia_xyz,'("He  ",3f18.13)') x,y,z
         end if
           if (atoms_to_remove_type(is).eq.71) then
            write(f_uvsia_xyz,'("H  ",3f18.13)') x,y,z
            write(f_uvbsia_xyz,'("H  ",3f18.13)') x,y,z
            write(f_uvbcsia_xyz,'("H  ",3f18.13)') x,y,z
         end if
         endif
      end do
        end do

         write(*,*) xmin,ymin,zmin,xmax,ymax,zmax

       ! then write the centers of C15 clusters ...
!       do i =1,nsia
!
!          write(*,*) i, rxb(nsite(i))/(a0/2.d0),ryb(nsite(i))/(a0/2.d0),rzb(nsite(i))/(a0/2.d0)
!
!       end do



       do is=1,nsia
         i=nsite(is)
         if (i.ne.0) then
           write(f_geom_org_data,'(a3,3(f18.13),"   1")')  '   ', rxb(i)*volc,ryb(i)*volc,rzb(i)*volc
           icnt=icnt+1
           write(f_sia_xyz,'("Cu  ",3f12.5)') rxb(i)*volc,ryb(i)*volc,rzb(i)*volc
           write(f_jsia_xyz,'("Cu  ",3f12.5)') rxb(i)*volc,ryb(i)*volc,rzb(i)*volc
           write(f_unit_cell,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
           write(f_formd_new_str,'(a3,3(f18.13))') 'Fe ', rxb(i)*volc,ryb(i)*volc,rzb(i)*volc
           write(f_geom_new_data,'(a3,3(f18.13))') '   ', rxb(i)*volc,ryb(i)*volc,rzb(i)*volc
           xp(1:3)=(/ rxb(i),ryb(i),rzb(i) /)*volc
           xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
           write(f_new_ndm_scl,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
         end if
       end do




      i_count=0
      do i=gndatoms+1,new_solid

        dist_ini=99999999.0
        do is=1,nsia
         it_a=nsite(is)
         if (it_a==0) it_a=i
         dist_current=dsqrt((rx(i)-rxb(it_a))**2+(ry(i)-ryb(it_a))**2+(rz(i)-rzb(it_a))**2)
           if (dist_current < dist_ini) dist_ini=dist_current

        end do
        if (dist_ini  > 0.1) then
           !
            i_count=i_count+1
            !write(*,*) rx(i), ry(i), rz(i), volc
            write(f_geom_org_data,'(a3,3(f18.13),"   1")')  '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
            icnt=icnt+1
            write(f_sia_xyz,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc
            write(f_jsia_xyz,'("Cu  ",3f12.5)') rx(i)*volc,ry(i)*volc,rz(i)*volc
            write(f_unit_cell,10) icnt,rx(i)*volc,ry(i)*volc,rz(i)*volc,in,in,in
            write(f_formd_new_str,'(a3,3(f18.13))') 'Fe ', rx(i)*volc,ry(i)*volc,rz(i)*volc
            write(f_geom_new_data,'(a3,3(f18.13))') '   ', rx(i)*volc,ry(i)*volc,rz(i)*volc
             xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
             xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
            write(f_new_ndm_scl,'(a3,3(f18.13),"   1")') '   ', xc(1:3)
           !
        end if

      end do

!                write(*,*) gndatoms, nsia, i_count
!            stop




       icnt=0
       do i=1,ntemp
        icnt =icnt+1
        x = rxb(i)*volc
        y = ryb(i)*volc
        z = rzb(i)*volc
        if (x.gt.xmin) then
         if (y.gt.ymin) then
          if (z.gt.zmin) then
           if (x.lt.xmax) then
            if (y.lt.ymax) then
             if (z.lt.zmax) then
                write(f_uvbcsia_xyz,'("Cu  ",3f18.13)') x,y,z
                lremove=.false.
                do is=1,ndatoms
                 if (i==nint(atoms_to_remove(is)))  lremove=.true.
              end do
                 if (.not.lremove) write(f_uvbcsia_sim_xyz,'("Cu  ",3f18.13)') x,y,z



             end if
            end if
           end if
          end if
         end if
        end if
       end do


       open(14,file='unit_vect',form='formatted')
       write (14,*) ' number of unit vectors 0 1 2 ou 3'
       write(14,*) 3
       write(14,*) 'unit vectors:'
       write(14,'(3f18.13)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
       write(14,'(3f18.13)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
       write(14,'(3f18.13)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

       write(f_usia_xyz,*) ' '
       write(f_usia_xyz,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
       write(f_usia_xyz,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
       write(f_usia_xyz,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc

       write(f_sia_xyz,*) ' '
       write(f_sia_xyz,'(3f12.5)') a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
       write(f_sia_xyz,'(3f12.5)') a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
       write(f_sia_xyz,'(3f12.5)') a(3,1)*volc, a(3,2)*volc, a(3,3)*volc




       close(f_sia_xyz)
       close(f_usia_xyz)
       close(f_jsia_xyz)
       close(f_unit_cell)
       close(14)



       return
       end



    subroutine write_formd_new_org (ntemp,nsize,rx,ry,rz,a,a0,volc)
    use math
    use pbc
    implicit none
    integer, intent(in)  :: ntemp,nsize
    real(8), intent(in)  :: rx(nsize),ry(nsize),rz(nsize)
    real(8), intent(in)  ::a(3,3),a0,volc
    real(8)   :: a02
    integer ::  i
    real(kind(1.d0)), dimension(3) :: xp,xc


              a02=a0/2.d0
 open (21,file='formd_new_ndm_org.str',form='formatted',status='unknown')

        write(21,'("1   1   1")')
      write(21,'(3f18.13)')a(1,1)*volc, a(1,2)*volc, a(1,3)*volc
        write(21,'(3f18.13)')a(2,1)*volc, a(2,2)*volc, a(2,3)*volc
        write(21,'(3f18.13)')a(3,1)*volc, a(3,2)*volc, a(3,3)*volc
      write(21,'(i8)') ntemp

       do i=1,ntemp
        xp(1:3)=(/ rx(i),ry(i),rz(i) /)*volc
        xc(1:3)=MatMul(inv_at(:,:)/volc,xp(:))
       write(21,'("   ",3(f18.13),"   1")') xc(1:3)
       end do

      close(21)
    return
    end





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
             write(*,*) '........', is_tot


return
end subroutine solid_and_defect_cut
