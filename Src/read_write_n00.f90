    integer ::      nat,ixstep,iystep,izstep
    integer ::   n1,n2,n3,i
    integer, allocatable, dimension(:) :: n1s, n2s,n3s
    
    open(unit=10,file='no00.inp')
    open(unit=11,file='no00_new.inp')

    read(10,*) nat
    allocate (n1s(nat),n2s(nat),n3s(nat))

     do i=1,nat
      read(10,  *) n1s(i),n2s(i),n3s(i)
     end do  

     write(*,*) '1 MAX MIN  AVERAGE', MAXVAL(n1s), MINVAL(n1s), int( (MAXVAL(n1s)+MINVAL(n1s))/2) 
     write(*,*) '2 MAX MIN  AVERAGE', MAXVAL(n2s), MINVAL(n2s), int( (MAXVAL(n2s)+MINVAL(n2s))/2) 
     write(*,*) '3 MAX MIN  AVERAGE', MAXVAL(n3s), MINVAL(n3s), int( (MAXVAL(n3s)+MINVAL(n3s))/2) 
     rewind(10)

    
    write(*,*) 'Type of no00 file ? (1 - traditional loops, 2 - C15 clusters, 3 - migration with care)' 
    read(*,*) itype
    if (itype==3)  then
     write(*,*) 'Not implemented yet'
     stop
    end if
 
    write(*,*) 'The step on X?'
    read(*,*) ixstep
    write(*,*) 'The step on Y?'
    read(*,*) iystep
    write(*,*) 'The step on Z?'
    read(*,*) izstep

    read(10,*) nat
    write(11,'(i7)') nat    
     do i=1,nat
      read(10,  *) n1,n2,n3
      write(11,'(3i7)') n1+ixstep,n2+iystep,n3+izstep
     end do  
    
     if (itype==2) then
     do i=1,nat
     read(10,*) n1
     write(11,'(i3)') n1
     end do


     end if 
    
    
    
end
