  subroutine read_tsia (tsia)
  integer :: tsia

       write(*,*) 'Choice the sia (1 - 110)(2 - 111)(21 - 111 with care)(22 - 111 m)(221 - 111 m with care)(3 - 100) (31 - 100 with care) (4-gao)(5 -vac)(6-octa) (67 - octa, twisted):'
       read(*,*) tsia
       write(*,*) 'excelent choice master'
       
       if (.not.((tsia.lt.7).or.(tsia==22).or.(tsia==21).or.(tsia==31).or.(tsia==67).or.(tsia==221))) then
        write(*,*) "This software it isi not designed for tsia", tsia
        write(*,*) "Try.....: 1, 2 , 21, 22, 221, 3, 31, 4, 5 , 6" 
        write(*,*) "Stop ...: read_tsia"
          stop
       end if 
  return

  end subroutine read_tsia