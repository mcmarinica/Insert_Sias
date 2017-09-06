 module migration_sia
   integer :: nfix,nmig
   integer, allocatable,dimension(:) :: type_migration,nsite_fixed,  &
                                        nsite_removed,nsite_trans,   &
                                        nsite_migra

   type migsia
   real(8)  :: atom(3,3)  ! atom(xyz,number) 
   end type migsia
 
   type(migsia), allocatable, dimension(:) :: debatoms,finatoms 
contains 


 subroutine allocate_migration (nsia)
   integer :: nsia


  allocate (type_migration(nsia))
  allocate (debatoms(nmig),finatoms(nmig),           &
            nsite_fixed(nfix),nsite_removed(2*nmig),  &
            nsite_migra(nmig),nsite_trans(nmig))
      




  end subroutine allocate_migration



end module migration_sia
