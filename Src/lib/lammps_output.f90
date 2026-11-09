!=======================================================================
! lammps_output.f90
!
! LAMMPS data files: atom style "atomic", one atom type, lengths in
! Angstrom (units metal). The files are read with the LAMMPS command
! read_data.
!
!   lammps_write_all   writes a complete file from arrays of coordinates
!   lammps_begin       opens a file and writes its header
!   lammps_atom        appends one atom
!   lammps_end         closes the file and checks the number of atoms
!
! The cell vectors are the rows of the matrix "cell", as in formd.str.
! LAMMPS needs the first vector along x and the second in the xy plane.
! A cell that does not satisfy this is not written: a warning is printed
! and the three routines that follow lammps_begin do nothing.
!=======================================================================
module lammps_output

  implicit none
  private
  public :: lammps_write_all, lammps_begin, lammps_atom, lammps_end

  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: tol = 1.d-10     ! a cell component below this is zero

  integer :: lmp_unit = -1
  integer :: natoms_announced = 0
  integer :: natoms_written = 0
  logical :: lmp_open = .false.
  character(len=256) :: lmp_name = ' '

contains

!-----------------------------------------------------------------------
  subroutine lammps_begin (filename, natoms, cell, title)
!   Opens "filename" and writes the header: title, number of atoms,
!   box. The atoms follow, one call of lammps_atom each.

    character(len=*), intent(in) :: filename   ! name of the data file
    integer, intent(in)          :: natoms     ! number of atoms that will be written
    real(dp), intent(in)         :: cell(3,3)  ! cell vectors, one per row
    character(len=*), intent(in) :: title      ! text of the first line

    logical :: orthogonal

    if (lmp_open) then
      write(*,'(4a)') 'ERROR <lammps_begin>: ', trim(lmp_name), &
                      ' is still open when starting ', trim(filename)
      stop
    end if

    if ( (abs(cell(1,2)) > tol) .or. (abs(cell(1,3)) > tol) .or.  &
         (abs(cell(2,3)) > tol) .or. (cell(1,1) <= tol) .or.      &
         (cell(2,2) <= tol) .or. (cell(3,3) <= tol) ) then
      write(*,'(3a)') 'WARNING <lammps_begin>: ', trim(filename), ' is not written.'
      write(*,'(a)')  '   LAMMPS needs the first cell vector along x and the second in the xy plane.'
      return
    end if

    orthogonal = (abs(cell(2,1)) <= tol) .and. (abs(cell(3,1)) <= tol) .and. &
                 (abs(cell(3,2)) <= tol)

    open (newunit=lmp_unit, file=filename, form='formatted', status='replace')
    lmp_open = .true.
    lmp_name = filename
    natoms_announced = natoms
    natoms_written = 0

    write(lmp_unit,'(2a)') '# ', trim(title)
    write(lmp_unit,'(a)')
    write(lmp_unit,'(i0,a)') natoms, ' atoms'
    write(lmp_unit,'(a)') '1 atom types'
    write(lmp_unit,'(a)')
    write(lmp_unit,'(2f24.16,a)') 0.d0, cell(1,1), '   xlo xhi'
    write(lmp_unit,'(2f24.16,a)') 0.d0, cell(2,2), '   ylo yhi'
    write(lmp_unit,'(2f24.16,a)') 0.d0, cell(3,3), '   zlo zhi'
    if (.not. orthogonal) then
      write(lmp_unit,'(3f24.16,a)') cell(2,1), cell(3,1), cell(3,2), '   xy xz yz'
    end if
    write(lmp_unit,'(a)')
    write(lmp_unit,'(a)') 'Atoms # atomic'
    write(lmp_unit,'(a)')

  end subroutine lammps_begin

!-----------------------------------------------------------------------
  subroutine lammps_atom (x, y, z)
!   Appends one atom, of type 1, at the Cartesian position (x,y,z).

    real(dp), intent(in) :: x, y, z

    if (.not. lmp_open) return
    natoms_written = natoms_written + 1
    write(lmp_unit,'(i10,i6,3f24.16)') natoms_written, 1, x, y, z

  end subroutine lammps_atom

!-----------------------------------------------------------------------
  subroutine lammps_end ()
!   Closes the file and compares the number of atoms written with the
!   number announced in the header.

    if (.not. lmp_open) return
    close (lmp_unit)
    lmp_open = .false.

    if (natoms_written == natoms_announced) then
      write(*,'(3a,i0,a)') 'LAMMPS data file written ..........: ', trim(lmp_name), &
                           ' (', natoms_written, ' atoms)'
    else
      write(*,'(3a)')      'ERROR <lammps_end>: ', trim(lmp_name), ' is inconsistent.'
      write(*,'(a,i0,a,i0,a)') '   Its header announces ', natoms_announced, &
                               ' atoms and ', natoms_written, ' were written.'
    end if

  end subroutine lammps_end

!-----------------------------------------------------------------------
  subroutine lammps_write_all (filename, natoms, cell, x, y, z, title)
!   Writes a complete file for the atoms (x(i), y(i), z(i)), i=1,natoms.

    character(len=*), intent(in) :: filename
    integer, intent(in)          :: natoms
    real(dp), intent(in)         :: cell(3,3)
    real(dp), intent(in)         :: x(natoms), y(natoms), z(natoms)
    character(len=*), intent(in) :: title

    integer :: i

    call lammps_begin (filename, natoms, cell, title)
    do i = 1, natoms
      call lammps_atom (x(i), y(i), z(i))
    end do
    call lammps_end ()

  end subroutine lammps_write_all

end module lammps_output
