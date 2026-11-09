!=======================================================================
! sites.f90
!
! Lattice sites: integer coordinates of the atoms, in units of a0/2, and
! search of the atoms that sit on the listed sites.
!
!   write_n123                 integer coordinates of all atoms; periodicity of the cell
!   analyze_config_sia         neighbour statistics of the listed sites (fort.567)
!   check_site                 atom of each listed site, and the atoms that a Gao, C15 or icosahedral unit removes
!   check_site_migration       the same for dumbbell jumps
!   care_check_site_migration  row neighbours of the jumping dumbbells (type 221)
!   fix_ico_type               orientation of each icosahedral unit, from the parity of its centre
!   fix_ico_type_mixed         the same for type 77
!=======================================================================

   subroutine write_n123(ntemp,nsize,a0,n1,n2,n3,rx,ry,rz,  &
                         n1MAX,n2MAX,n3MAX)
   use math
   use the_111_smooth, ONLY : pbc_111
   use the_100_smooth, ONLY : pbc_100
   implicit none

   integer :: nsize, i,n1MAX,n2MAX,n3MAX,ntemp
   integer  :: n1(nsize),n2(nsize),n3(nsize)
   real(8)     :: rx(nsize),ry(nsize),rz(nsize)
   integer :: itemp,ngcd,n1pbc,n2pbc,n3pbc

   double precision     :: a0,a02

   a02= a0/2.d0

     do i=1,ntemp
          n1(i)=nint(rx(i)/a02)
          n2(i)=nint(ry(i)/a02)
          n3(i)=nint(rz(i)/a02)
     end do
    n1MAX = maxval(n1(1:ntemp))
    n2MAX = maxval(n2(1:ntemp))
    n3MAX = maxval(n3(1:ntemp))

    n1pbc=n1MAX+1
    n2pbc=n2MAX+1
    n3pbc=n3MAX+1

    print*, 'n1MAX,n2MAX,n3MAX',n1MAX,n2MAX,n3MAX

    call gcd (n1pbc,n2pbc,itemp)
    call gcd (itemp,n3pbc,ngcd)
    pbc_111=n1pbc*(n2pbc/ngcd)*(n3pbc/ngcd)
    pbc_100= n3pbc
     write(*,*) 'We will have periodicity in <111> direction in', pbc_111
     write(*,*) 'We will have periodicity in <100> direction in', pbc_100

   return
   end

!=======================================================================

   subroutine analyze_config_sia (tsia,nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)

   use defect_types
   implicit none
   integer  :: nsia, n1MAX,n2MAX,n3MAX,i,j
   integer  :: ns1(nsia),ns2(nsia),ns3(nsia)
   integer  :: ii,neigh(5),tsia
   real(8)  :: bcc(5),dij,delta=0.1d0

  if (.not.((tsia==DEFECT_SIA_110).or.(tsia==DEFECT_SIA_111).or. &
            (tsia==DEFECT_SIA_100).or.(tsia==DEFECT_SIA_111_CARE))) then
   write(*,*) 'WARNING ... the analyze_config SIA .... Not yet implemented for tsia  ',tsia
   !stop
  end if

  ! the neighbour shells are defined for types 2 and 3 only
  bcc(:)=0.d0

  if (tsia==DEFECT_SIA_111) then
 ! in a0 units the bcc NN distances
  bcc(1)=dsqrt(3.d0)/2.d0
  bcc(2)=1.d0
  bcc(3)=dsqrt(2.d0)
  bcc(4)=dsqrt(3.d0)
  bcc(5)=2.d0
 end if

if (tsia==DEFECT_SIA_100) then
  bcc(1)=dsqrt(3.d0)/2.d0
  bcc(2)=1.d0
  bcc(3)=dsqrt(2.d0)
  bcc(4)=dsqrt(11.d0)/2.d0
  bcc(5)=2.d0
end if

 ! in internal units of a0/2 this distances gives
  bcc(:)=2.d0*bcc(:)

 write(567,'(i8)')  nsia
 do i=1,nsia
    neigh(:) = 0
  do j=1,nsia
  ! if (tsia==DEFECT_SIA_100) then
  !  dij=dsqrt(dble(                      &
  !            (ns1(i)-ns1(j))**2 +  &
  !            (ns2(i)-ns2(j))**2    &
  !            ))
  !
  !else
    dij=dsqrt(dble(                      &
              (ns1(i)-ns1(j))**2 +  &
              (ns2(i)-ns2(j))**2 +  &
              (ns3(i)-ns3(j))**2    &
              ))
   !end if
    do ii=1,5
     if ((dij<=(bcc(ii)+delta)).and.(dij>=(bcc(ii)-delta))) then
        neigh(ii)=neigh(ii)+1
     end if
    end do
   end do !j-nsia
   write(567,'(i8,5i6)') i,neigh(1:5)
 end do  !i-nsia

   return
   end

!=======================================================================

       subroutine check_site (tsia,ntemp,nsize,nsia,nlac_remove,n1,n2,n3,ns1,ns2,ns3, &
                              nsite,nsite_m,nsite_octa_vac_remove,nsite_octa_vac_remove_type,type_octa)

       use octahedre
       use icosaedre
       use the_111_smooth
       use the_100_smooth
       use defect_types
       implicit none
       integer :: ntemp,itst,itst2,ii,jj,i,tsia,is,nsia,nsize,ti_octa,ilac_remove
       integer, intent(out)   :: nlac_remove
       integer :: minval1,minval2,minval3
       integer, dimension(nsize)  :: n1,n2,n3
       integer, dimension(nsia)  :: nsite,ns1,ns2,ns3,type_octa
       integer, dimension(4,nsia) :: nsite_m,ns1_m,ns2_m,ns3_m
       integer, dimension(10*nsia) :: nsite_octa_vac_remove, ns1_octa,ns2_octa,ns3_octa
       integer, dimension(10*nsia) :: nsite_octa_vac_remove_type
       integer :: ti_icos

        itst=0
        itst2=0
        do is=1,nsia
         !debug write(*,*) 'ns1',ns1(1:3)
          do i=1,ntemp
             if (n1(i)==ns1(is)) then
              if (n2(i)==ns2(is)) then
               if (n3(i)==ns3(is)) then
                itst =itst +1
                nsite(is)= i
               end if
              end if
             end if

          if (tsia==DEFECT_SIA_111_CARE) then
           do jj=-2,2,1
            !
             if (n1(i)==ns1_111(jj,is)) then
              if (n2(i)==ns2_111(jj,is)) then
               if (n3(i)==ns3_111(jj,is)) then
                itst2 =itst2 +1
                nsite_111(jj,is)= i
               end if
              end if
             end if
            !
            end do ! jj
          end if !tsia==DEFECT_SIA_111_CARE

          if (tsia==DEFECT_SIA_100_CARE) then
           do jj=-2,2,1
            !
             if (n1(i)==ns1_100(jj,is)) then
              if (n2(i)==ns2_100(jj,is)) then
               if (n3(i)==ns3_100(jj,is)) then
                itst2 =itst2 +1
                nsite_100(jj,is)= i
               end if
              end if
             end if
            !
            end do ! jj
          end if !tsia==DEFECT_SIA_100_CARE



        end do !ntemp
      end do   !is
      if (.not. defect_is_ico(tsia)) then
        if (itst.ne.nsia) then
         print*, 'there is a big problem in input'
         print*, 'Check check_site'
         print*, '< check_site >'
         print*, itst, nsia
         stop
       end if
      end if


       if (tsia==DEFECT_GAO) then

       do ii=1,nsia
          nsite_m(1,ii)=nsite(ii)

           ns1_m(1,ii) = ns1(ii)
           ns1_m(2,ii) = ns1(ii)+1
           ns1_m(3,ii) = ns1(ii)-1
           ns1_m(4,ii) = ns1(ii)+1

           ns2_m(1,ii) = ns2(ii)
           ns2_m(2,ii) = ns2(ii)-1
           ns2_m(3,ii) = ns2(ii)-1
           ns2_m(4,ii) = ns2(ii)+1

           ns3_m(1,ii) = ns3(ii)
           ns3_m(2,ii) = ns3(ii)-1
           ns3_m(3,ii) = ns3(ii)+1
           ns3_m(4,ii) = ns3(ii)+1
         end do


        itst=0
        do is=1,nsia
           do i=1,ntemp
            do ii=2,4
               if (n1(i)==ns1_m(ii,is)) then
                if (n2(i)==ns2_m(ii,is)) then
                 if (n3(i)==ns3_m(ii,is)) then
                   itst =itst +1
                   nsite_m(ii,is)= i
                 end if
                end if
               end if
             end do
            end do
          end do

        if (itst.ne.3*nsia) then
          print*, 'there is a big problem in input'
          print*, 'Check check_site 3*'
          print*, itst, nsia
          stop
          end if

        !C15 case ...
       else if (defect_is_c15(tsia)) then



       ilac_remove=0
       do ii=1,nsia
         ti_octa=type_octa(ii)
         if ( (ti_octa.eq.1).or.(ti_octa.eq.2) ) then  ! This is perfect C15  I2
           do jj=1,4
             ilac_remove=ilac_remove+1
             ns1_octa(ilac_remove) = ns1(ii)+nint(octo(ti_octa)%atom(1,jj))
             ns2_octa(ilac_remove) = ns2(ii)+nint(octo(ti_octa)%atom(2,jj))
             ns3_octa(ilac_remove) = ns3(ii)+nint(octo(ti_octa)%atom(3,jj))
             nsite_octa_vac_remove_type(ilac_remove)= 70  ! 70 = gao
             !debug write(*,'(2i5,3i3)') jj, ti_octa,&
             !debug nint(octo(ti_octa)%atom(1,jj)), nint(octo(ti_octa)%atom(2,jj)), nint(octo(ti_octa)%atom(3,jj))
           end do

           do jj=1,6
             ilac_remove=ilac_remove+1
             ns1_octa(ilac_remove) = ns1(ii)+nint(octo(ti_octa)%removeatom(1,jj))
             ns2_octa(ilac_remove) = ns2(ii)+nint(octo(ti_octa)%removeatom(2,jj))
             ns3_octa(ilac_remove) = ns3(ii)+nint(octo(ti_octa)%removeatom(3,jj))
             nsite_octa_vac_remove_type(ilac_remove)= 71  ! 71 = sia
             !debug write(*,'(2i5,3i3)') jj, ti_octa,&
             !debug  nint(octo(ti_octa)%removeatom(1,jj)), nint(octo(ti_octa)%removeatom(2,jj)), nint(octo(ti_octa)%removeatom(3,jj))
           end do

          else if (ti_octa.eq.4) then  ! This is the new object for odd clusters I2 + I
            ti_octa=ti_octa-2
            ilac_remove=ilac_remove+1
            ns1_octa(ilac_remove) = ns1(ii)+nint(octo(ti_octa)%atom(1,1))
            ns2_octa(ilac_remove) = ns2(ii)+nint(octo(ti_octa)%atom(2,1))
            ns3_octa(ilac_remove) = ns3(ii)+nint(octo(ti_octa)%atom(3,1))
            nsite_octa_vac_remove_type(ilac_remove)= 70  ! 70 = gao


            ilac_remove=ilac_remove+1

            ns1_octa(ilac_remove) = ns1(ii)+nint(octo(ti_octa)%removeatom(1,4))
            ns2_octa(ilac_remove) = ns2(ii)+nint(octo(ti_octa)%removeatom(2,4))
            ns3_octa(ilac_remove) = ns3(ii)+nint(octo(ti_octa)%removeatom(3,4))
            nsite_octa_vac_remove_type(ilac_remove)= 71  ! 71 = sia
          end if

        end do
        !ico case ...
       else if (tsia==DEFECT_ICO) then
          ilac_remove=0
          do ii=1,nsia
            ti_icos=1
              do jj=1,6
                ilac_remove=ilac_remove+1
                ns1_octa(ilac_remove) = ns1(ii)+nint(icos(ti_icos)%removeatom(1,jj))
                ns2_octa(ilac_remove) = ns2(ii)+nint(icos(ti_icos)%removeatom(2,jj))
                ns3_octa(ilac_remove) = ns3(ii)+nint(icos(ti_icos)%removeatom(3,jj))
                nsite_octa_vac_remove_type(ilac_remove)= 70  ! 70 = Sias on the face of the cube
                !debug write(*,'(2i5,3i3)') jj, ti_octa,&
                !debug nint(octo(ti_octa)%atom(1,jj)), nint(octo(ti_octa)%atom(2,jj)), nint(octo(ti_octa)%atom(3,jj))
              end do

              do jj=1,1
                ilac_remove=ilac_remove+1
                ns1_octa(ilac_remove) = ns1(ii)+nint(icos(ti_icos)%add(1,jj))
                ns2_octa(ilac_remove) = ns2(ii)+nint(icos(ti_icos)%add(2,jj))
                ns3_octa(ilac_remove) = ns3(ii)+nint(icos(ti_icos)%add(3,jj))
                nsite_octa_vac_remove_type(ilac_remove)= 71  ! 71 =  the center of the icosaedre
                !debug write(*,'(2i5,3i3)') jj, ti_octa,&
                !debug  nint(octo(ti_octa)%removeatom(1,jj)), nint(octo(ti_octa)%removeatom(2,jj)), nint(octo(ti_octa)%removeatom(3,jj))
             end do
           end do

       else if (tsia==DEFECT_ICO_DUMBBELLS) then
          ilac_remove=0
          do ii=1,nsia_ico
            ti_icos=1
              do jj=1,6
                ilac_remove=ilac_remove+1
                ns1_octa(ilac_remove) = ns1(ii)+nint(icos(ti_icos)%removeatom(1,jj))
                ns2_octa(ilac_remove) = ns2(ii)+nint(icos(ti_icos)%removeatom(2,jj))
                ns3_octa(ilac_remove) = ns3(ii)+nint(icos(ti_icos)%removeatom(3,jj))
                nsite_octa_vac_remove_type(ilac_remove)= 70  ! 70 = Sias on the face of the cube
                !debug write(*,'(2i5,3i3)') jj, ti_octa,&
                !debug nint(octo(ti_octa)%atom(1,jj)), nint(octo(ti_octa)%atom(2,jj)), nint(octo(ti_octa)%atom(3,jj))
              end do

              do jj=1,1
                ilac_remove=ilac_remove+1
                ns1_octa(ilac_remove) = ns1(ii)+nint(icos(ti_icos)%add(1,jj))
                ns2_octa(ilac_remove) = ns2(ii)+nint(icos(ti_icos)%add(2,jj))
                ns3_octa(ilac_remove) = ns3(ii)+nint(icos(ti_icos)%add(3,jj))
                nsite_octa_vac_remove_type(ilac_remove)= 71  ! 71 =  the center of the icosaedre
                !debug write(*,'(2i5,3i3)') jj, ti_octa,&
                !debug  nint(octo(ti_octa)%removeatom(1,jj)), nint(octo(ti_octa)%removeatom(2,jj)), nint(octo(ti_octa)%removeatom(3,jj))
             end do
           end do
        end if !if tsia eq 4

       if (defect_is_cluster(tsia)) then

         nlac_remove=ilac_remove
         minval1=MINVAL(ns1_octa(1:nlac_remove))
         minval2=MINVAL(ns2_octa(1:nlac_remove))
         minval3=MINVAL(ns3_octa(1:nlac_remove))
         write(*,'("The minval on x,y and z.....:",3i6)') &
                    minval1,minval2,minval3
         if(  (minval1.lt.0).or.(minval2.lt.0).or.(minval3.lt.0) ) then
           Write(*,*) &
            'There are atoms which exceed the box .... change the input'
          stop
         end if

         itst=0
         do is=1,nlac_remove
           do i=1,ntemp
             if (n1(i)==ns1_octa(is)) then
               if (n2(i)==ns2_octa(is)) then
                 if (n3(i)==ns3_octa(is)) then
                   itst =itst +1
                   nsite_octa_vac_remove(itst)= i   ! This site will be removed in the end
                                                  ! because it is a Gao or a vacancy in a dumbbell.
                   !debug write(*,*) 'check',is, itst, nsite_octa_vac_remove(itst)
                 end if
               end if
             end if
           end do
         end do
         if (.not. defect_is_ico(tsia)) then
         if (itst.ne.ilac_remove) then
           print*, 'there is a big problem in input'
           print*, 'Check <check_site> 10*'
           print*, itst, nsia
           !stop
         end if
         end if

         nlac_remove=itst
       end if   !the last if
       return
       end

!=======================================================================

       subroutine check_site_migration (ntemp,nsize,nsia,       &
                                        n1,n2,n3,                    &
                                        ns1,ns2,ns3)

       use migration_sia
       implicit none
       integer :: ntemp,itst,i,is,nsia,nsize
       integer, dimension(nsia)  :: nsite,ns1,ns2,ns3
       integer, dimension(nsize)  :: n1,n2,n3
       integer   :: ifixed,iremoved,imigra1,imigra2

! check time ............

        itst=0
        do is=1,nsia
         !debug write(*,*) 'ns1',ns1(1:3)
          do i=1,ntemp
             if (n1(i)==ns1(is)) then
              if (n2(i)==ns2(is)) then
               if (n3(i)==ns3(is)) then
                itst =itst +1
                nsite(is)= i
               end if
              end if
             end if
          end do
         end do

        if (itst.ne.nsia) then
         print*, 'there is a big problem in input'
         print*, 'Check check_site'
         print*, '< check_site >'
         print*, itst, nsia
         stop
       end if


! all the atoms which should be removed ...
        ifixed=0        ! a counter for non-migrating atoms (or fixed)
        imigra1=0
        imigra2=0
        iremoved=0      ! the atoms which should be removed from the initial and final file.
        do is=1,nsia
          do i=1,ntemp

            if (n1(i)==ns1(is)) then
              if (n2(i)==ns2(is)) then
               if (n3(i)==ns3(is)) then
                 !
                 if (type_migration(is)==1) then
                  iremoved=iremoved+1
                  imigra1=imigra1+1
                  nsite_removed (iremoved)= i
                  nsite_migra (imigra1)= i
                 end if
                 if (type_migration(is)==0) then
                   ifixed=ifixed+1
                   nsite_fixed  (ifixed )= i
                 end if
                 !

              end if
           end if
          end if

           if (type_migration(is)==1) then
             if (n1(i)==ns1(is)+1) then
              if (n2(i)==ns2(is)+1) then
               if (n3(i)==ns3(is)+1) then
                iremoved=iremoved+1
                imigra2=imigra2+1
                nsite_removed (iremoved)= i
                nsite_trans (imigra2) = i
               end if
              end if
             end if
           end if

          end do   !ntemp
       end do      !nsia

       if (nfix/=ifixed) then
          write(*,*) "Problem in nfix and ifixed",nfix,ifixed
          stop
       end if

       if (nmig /= iremoved/2) then
        write(*,*) "Problem in nmig and iremoved", nmig,iremoved
       end if



       return
       end

!=======================================================================

  subroutine care_check_site_migration (ntemp,nsize,nsia,       &
                                        n1,n2,n3)

  use migration_sia
  use the_111_smooth, ONLY: deb_ns1_111,deb_ns2_111,deb_ns3_111,   &
                                 fin_ns1_111,fin_ns2_111,fin_ns3_111,   &
                            deb_nsite_111,fin_nsite_111
  implicit none
  integer, intent (in) :: ntemp,nsia,nsize
  integer :: i,is,jj
  integer, dimension(nsize), intent(in)  :: n1,n2,n3
  integer   :: itst1,itst2

! check time ............



! all the atoms which should be removed ...
   itst1=0        ! a counter for non-migrating atoms (or fixed)
   itst2=0
   do is=1,nsia
     do i=1,ntemp
      do jj=-2,2,1
       !
        if (n1(i)==deb_ns1_111(jj,is)) then
         if (n2(i)==deb_ns2_111(jj,is)) then
          if (n3(i)==deb_ns3_111(jj,is)) then
           itst2 =itst2 +1
           deb_nsite_111(jj,is)= i
          end if
         end if
        end if
       !
       end do ! jj

      do jj=-2,2,1
       !
        if (n1(i)==fin_ns1_111(jj,is)) then
         if (n2(i)==fin_ns2_111(jj,is)) then
          if (n3(i)==fin_ns3_111(jj,is)) then
           itst2 =itst2 +1
            fin_nsite_111(jj,is)= i
          end if
         end if
        end if
       !
       end do ! jj

     end do   !ntemp
  end do      !nsia

  !if (nfix/=ifixed) then
  !   write(*,*) "Problem in nfix and ifixed",nfix,ifixed
  !   stop
  !end if

  !if (nmig /= iremoved/2) then
  ! write(*,*) "Problem in nmig and iremoved", nmig,iremoved
  !end if



  return
end subroutine   care_check_site_migration

!=======================================================================

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

!=======================================================================

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

!=======================================================================
