!=======================================================================
! input.f90
!
! Reading of the input: the defect type from standard input, the reference
! crystal (formd.str) and the list of defect sites (no00.inp, or no00m.inp
! for dumbbell jumps).
!
!   read_tsia                  asks for the defect type
!   read_ntemp                 opens formd.str and reads the number of atoms
!   read_nsia                  opens the list of sites and reads its first line
!   read_input                 reads the cell and the atoms of formd.str
!   read_config_sia            reads the sites
!   read_config_sia_migration  reads the sites and jump flags of no00m.inp
!   read_config_sia_octa       reads the types of the C15 units
!   read_config_sia_ico        reads the extra dumbbells and atoms of type 77
!=======================================================================

  subroutine read_tsia (tsia)
  use defect_types
  implicit none
  integer :: tsia

       write(*,*) 'Choice the sia (1 - 110)(2 - 111)(21 - 111 with care) (210-Frank loop) &
                  (22 - 111 m)(221 - 111 m with care)(3 - 100) (31 - 100 with care) &
                  (4-gao)(5 -vac)(6-octa) (67 - octa, twisted)(7 -ico)(77-ico, dumbells):'
       read(*,*) tsia
       write(*,*) 'excelent choice master'

       if (.not. defect_is_known(tsia)) then
        write(*,*) "This software it isi not designed for tsia", tsia
        write(*,*) "Try.....: 1, 2, 21, 210, 22, 221, 3, 31, 4, 5, 6, 67, 7, 77"
        write(*,*) "Stop ...: read_tsia"
          stop
       end if
  return

  end subroutine read_tsia

!=======================================================================

   subroutine read_ntemp (ntemp)
    integer :: ntemp

   open (10,file='formd.str',form='formatted',status='unknown')
     read (10,'(i8)') ntemp

    return
    end

!=======================================================================

   subroutine read_nsia (nsia,tsia)
   use icosaedre, only: nsia_ico, nsia_dumbell
   use migration_sia
   use defect_types
   implicit none
   integer :: nsia,tsia

   if (defect_is_migration(tsia)) then
    open (11,file='no00m.inp',form='formatted',status='unknown')
    read (11,*) nsia,nmig,nfix
    if ((nfix+nmig)/=nsia) then
      write(*,*) "nfix + nmig should be equal to nsia, nsia,nfix,nmig.....:",nsia,nfix,nmig
      write(*,*) "stop in read_nsia"
      stop
    end if

   else

    open (11,file='no00.inp',form='formatted',status='unknown')
    if (tsia==DEFECT_ICO_DUMBBELLS) then
       read (11,*) nsia_ico, nsia_dumbell
       nsia = nsia_ico
     else
       read (11,*) nsia
    end if
   end if

   return
   end

!=======================================================================

   subroutine read_input(ntemp,nsize,a,rx,ry,rz)
   use math
   use pbc
   implicit none
   integer              :: ntemp,nsize,i
   double precision     :: rx(nsize),ry(nsize),rz(nsize)
   double precision,dimension(:,:) :: a(3,3)
   character*3 :: sym


    read(10,'(3f18.13)')a(1,1), a(1,2), a(1,3)
    read(10,'(3f18.13)')a(2,1), a(2,2), a(2,3)
    read(10,'(3f18.13)')a(3,1), a(3,2), a(3,3)

    at(:,1)=a(1,:)
    at(:,2)=a(2,:)
    at(:,3)=a(3,:)
    ! Inverse matrix
    CALL MatInv(at,inv_at)

  do i=1,ntemp
        read(10,'(a3,3(f18.13))') sym, rx(i),ry(i),rz(i)
    end do



   return
   end

!=======================================================================

   subroutine read_config_sia (nsia,ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
   integer  :: nsia, n1MAX,n2MAX,n3MAX,i
   integer  :: ns1(nsia),ns2(nsia),ns3(nsia)


      do i=1,nsia
       read (11,*) ns1(i),ns2(i),ns3(i)
       !write(*,*) 'ns11',ns1(i),ns2(i),ns3(i)
        if ( (ns1(i) > n1MAX+1) .or.(ns2(i) > n2MAX+1) .or. &
          (ns3(i) > n3MAX+1) ) then
        print*, 'Big problem in read_config_sia'
        print*, 'The cell is not big enough. Increase the size!'
        print*, '<read_config_sia>', n1MAX+1, n2MAX+1, n3MAX+1
        stop
        end if
      end do
   return
   end

!=======================================================================

   subroutine read_config_sia_migration (nsia, ns1,ns2,ns3,n1MAX,n2MAX,n3MAX)
   use migration_sia
   integer  :: nsia, n1MAX,n2MAX,n3MAX,i
   integer  :: ns1(nsia),ns2(nsia),ns3(nsia)
   integer  :: i1,i2
     i1=0
     i2=0

      do i=1,nsia
       read (11,*) ns1(i),ns2(i),ns3(i),type_migration(i)
       !write(*,*) 'ns11',ns1(i),ns2(i),ns3(i)
        if ( (ns1(i) > n1MAX+1) .or.(ns2(i) > n2MAX+1) .or. &
          (ns3(i) > n3MAX+1) ) then
        print*, 'Big problem in read_config_sia'
        print*, 'The cell is not big enough. Increase the size!'
        print*, '<read_config_sia_migration>'
        stop
        end if
        if (.not.((type_migration(i)==0).or.(type_migration(i)==1))) then
          write(*,*) "The 4th colomn of no00m.inp should be 0 or 1"
          write(*,*) "stop in read_config_sia_migration"
          stop
        end if

        if (type_migration(i)==1) then
         if ( (ns1(i) > n1MAX-1) .or.(ns2(i) > n2MAX-1) .or. &
           (ns3(i) > n3MAX-1) ) then
           print*, 'Big problem in read_config_sia_migration'
           print*, 'The cell is not big enough. Increase the size!'
           print*, '<read_config_sia_migration>'
           stop
         end if
         if ( (ns1(i)==0) .or.(ns2(i)==0) .or. &
           (ns3(i)==0 )) then
           print*, 'Big problem in read_config_sia_migrration'
           print*, 'You cannot have 0 in the coordinates for migrating interstitials'
           print*, '<read_config_sia_migration>'
           stop
         end if
        end if


       if (type_migration(i)==1) i1=i1+1
       if (type_migration(i)==0) i2=i2+1
      end do

     if (i1/=nmig) then
      write(*,*) "The number of migrating atoms diff from nmig ...i1 and nmig ...:", i1,nmig
      stop
     end if

     if (i2/=nfix) then
      write(*,*) "The number of fixed atoms diff from nfix ...i2 and nfix ...:", i2,nfix
      stop
     end if



   return
   end

!=======================================================================

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

!=======================================================================

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

!=======================================================================
