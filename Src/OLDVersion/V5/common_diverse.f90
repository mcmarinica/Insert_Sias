!#####################################################################
       	  
       	  
       subroutine check_site (tsia,ntemp,nsize,nsia,nlac_remove,n1,n2,n3,ns1,ns2,ns3, &
                              nsite,nsite_m,nsite_octa_vac_remove,nsite_octa_vac_remove_type,type_octa)		         
       
       use octahedre 
       implicit none
       integer :: ntemp,itst,ii,jj,i,tsia,is,nsia,nsize,ti_octa,ilac_remove
       integer, intent(out)   :: nlac_remove
       integer :: minval1,minval2,minval3
       integer, dimension(nsize)  :: n1,n2,n3
       integer, dimension(nsia)  :: nsite,ns1,ns2,ns3,type_octa
       integer, dimension(4,nsia) :: nsite_m,ns1_m,ns2_m,ns3_m		
       integer, dimension(10*nsia) :: nsite_octa_vac_remove, ns1_octa,ns2_octa,ns3_octa		
       integer, dimension(10*nsia) :: nsite_octa_vac_remove_type

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


       if (tsia.eq.4) then	
        
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

        
        else if (tsia.eq.6) then

       	
       	
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
		   nsite_octa_vac_remove(is)= i   ! This site will be removed in the end
                                                  ! because it is a Gao or a vacancy in a dumbbell.
		   !debug write(*,*) 'check',i,n3(i),ii,is
		 end if 
		end if
	       end if
	       



	       
	    end do
	  end do  
	     			
        
	if (itst.ne.ilac_remove) then
	  print*, 'there is a big problem in input'
	  print*, 'Check <check_site> 10*'
	  print*, itst, nsia
	  stop	   	
  	end if	

	
	end if   !if tsia.eq.4
				
	return
	end			       

       	  
!#####################################################################

       subroutine neighb_sia (nall,nsia,volc,nsite,rx,  &
                       ry,rz,a)
       implicit none
       integer,   parameter :: NVMAX=6
       integer  :: nall,nsia,i,nsia2,j
       double precision     :: rx(nall),ry(nall),rz(nall)
       double precision     :: rs(3,nsia),ds(nsia*nsia),dsn(nsia*nsia)
       integer  :: nsite(nsia),indx(nsia*nsia)       
       double precision     ::  volc        
       double precision,dimension(:,:) :: a(3,3)
        integer          :: nvv(NVMAX), nvvf(NVMAX)
	double precision :: dv(NVMAX),dvf(NVMAX),bcc(NVMAX-1)        
      real(8)  :: bn
	 

        bcc(1)=1.d0
	bcc(2)=2.d0/dsqrt(3.d0)
	bcc(3)=dsqrt(2.d0)*2.d0/dsqrt(3.d0)
	bcc(4)=dsqrt(11.d0)/2.d0*bcc(2)
	bcc(5)=2.d0
	
	a(:,:)=a(:,:)*volc
	
        do i=1,nsia
	 rs(1,i)=rx(nsite(i))*volc
	 rs(2,i)=ry(nsite(i))*volc
	 rs(3,i)=rz(nsite(i))*volc
	end do

         
       call d_pbc(a,rs,ds,nsia)
       nsia2=nsia*nsia
       call indexx(nsia2,ds,indx) 
       do i=1,nsia2
        j=indx(i)
        dsn(i)=ds(j)
       end do 
       
       call countn(nsia2,NVMAX,dsn,nvv,dv)
       
      
      nvvf(:)=0
      nvvf(1)=nvv(1)
      do i=1,NVMAX-1
       bn=bcc(i)*dv(2)
       !debug_write_bn write(*,*) bn
       do j = 2,NVMAX
         if (dabs(dv(j)-bn).lt.0.001) then
	  nvvf(i+1)=nvv(j)
	  dvf(i+1)=dv(j)
	 end if
       end do
      end do      
       
 	do i=1,NVMAX
 	 write(*,'(2i6,f12.5)') i,nvv(i),dv(i)
 	end do

       write(*,'("=========================================")')
       write(*,'("      n1","      n2","      n3","      n4","      n5","     ntot")')
       write(*,'("=========================================")')
       write(*,'(6i8)')   (nvvf(i)/2,i=2,6),nvvf(1) 
       write(*,'(6f8.3)') (dvf(i),i=2,6   ),dvf(1)
       return
       end
!#####################################################################
       
       
      subroutine countn(NRMAX,NVMAX,a,nvv,dv)
! input the vector a(NRMAX)
! the output: the different values in dv(NVMAX)
!             how many values are for dv(i) recorded in nvv(i)     
      implicit none
      integer :: NRMAX,NVMAX
      double precision  :: at,ac,  a(NRMAX),dv(NVMAX)
      integer ::  nvv(NVMAX),ij,i,iv
      	     
      at=a(1)
      ij=1
      iv=1
        do i=2,NRMAX
	  ac=dabs(a(i)-at)
	  if (ac.lt.0.01) then
	   ij=ij+1
	  else
	   iv=iv+1   
            if (iv.gt.NVMAX) then
	     go to 113
	    end if 
	   ij=1
	  end if
	   nvv(iv)=ij
	   dv(iv)=a(i)

  	 !debug  write(*,'(2f8.2,2i6,f8.3)') a(i),ac,ij,iv,dv(iv)
	  at=a(i)
	 end do
	 
113     continue 
       return
       end     


!*********************************************************************!
        SUBROUTINE  d_pbc (unit_vect,tau,ds,nat_up)
        implicit none
        integer  nat_up,i,j,ni
	double precision x,y,z,a_side,b_side,c_side,a_side2,&
	                 b_side2,c_side2,a,b,c
        double precision tau(3,nat_up),     &
                          unit_vect(3,3),  &
			  dc(3,3),dci(3,3), &
			  ds( nat_up*nat_up )
!
!       determination du nombre de voisins jusqu'a Rmax
!	
	
!..............pbc setup

	a_side=sqrt(DOT_PRODUCT(unit_vect(:,1),unit_vect(:,1)))
	b_side=sqrt(DOT_PRODUCT(unit_vect(:,2),unit_vect(:,2)))
	c_side=sqrt(DOT_PRODUCT(unit_vect(:,3),unit_vect(:,3)))
	
	
	a_side2=a_side/2.d0
	b_side2=b_side/2.d0
	c_side2=c_side/2.d0
	
	dc(:,1)=unit_vect(:,1)/a_side
	dc(:,2)=unit_vect(:,2)/b_side
	dc(:,3)=unit_vect(:,3)/c_side
       
        call inverse(dc,dci)   

        ni=0
	Do  i=1,nat_up
         Do j=1,nat_up
            ni=ni+1
	    x = tau(1,i)-tau(1,j)          
            y = tau(2,i)-tau(2,j)          
            z = tau(3,i)-tau(3,j)   
             

            call fromto(x,y,z,a,b,c,1,dc,dci)

               if (a.lt.-a_side2    ) a=a+a_side
               if (a.gt. a_side2    ) a=a-a_side

               if (b.lt.-b_side2    ) b=b+b_side
               if (b.gt. b_side2    ) b=b-b_side

               if (c.lt.-c_side2    ) c=c+c_side
               if (c.gt. c_side2    ) c=c-c_side

            call fromto(x,y,z,a,b,c,-1,dc,dci)
    
           ds(ni)=dsqrt(x**2+y**2+z**2)


        end do
       end do
	   

       RETURN
       END
!*********************************************************************!


      subroutine fromto(x,y,z,a,b,c,imode,dc,dci)
!**********************************************************************
!*								      *
!* NAME: "FROMTO"						      *
!*								      *
!* USAGE: Transformation between crystallographic and Cartesian       *
!*	  coordinate frames.					      *
!*								      *
!* METHOD:			-1				      *
!*		[a]    [ax bx cx]  [x]       [x]		      *
!*		[b] =  [ay by cy]  [y] = DCI [y]		      *
!*		[c]    [az bz cz]  [z]       [z]		      *
!*								      *
!*		[x]    [ax bx cx]  [a]       [a]		      *
!*		[y] =  [ay by cy]  [b] = DC  [b]		      *
!*		[z]    [az bz cz]  [c]       [c]		      *
!*								      *
!* where matrix DC is a matrix of direction cosines of  	      *
!* crystallographic axes.					      *
!*								      *
!* INPUT/OUTPUT: 'x, y, z' - Cartesian coordinates;		      *
!*		 'a, b, c' - crystallographic coordinates;	      *
!*               'imode'   - direction of transformation:             *
!*                imode= 1  Cartesian --> crystallographic;           *
!*                imode=-1  crystallographic --> Cartesian.           *
!*                                                                    *                                                                    *
!* NOTE: Should be used together with (and after) MDCELL and MDCELL2. *
!**********************************************************************
      implicit none
      double precision :: dc(3,3),dci(3,3)
      integer :: imode
      real(8) :: a,b,c,x,y,z

      if (imode.eq.1) then

         a=dci(1,1)*x + dci(1,2)*y + dci(1,3)*z
         b=dci(2,1)*x + dci(2,2)*y + dci(2,3)*z
         c=dci(3,1)*x + dci(3,2)*y + dci(3,3)*z

      else if (imode.eq.-1) then

         x=dc(1,1)*a + dc(1,2)*b + dc(1,3)*c
         y=dc(2,1)*a + dc(2,2)*b + dc(2,3)*c
         z=dc(3,1)*a + dc(3,2)*b + dc(3,3)*c

      end if

      return
      end


!**********************************************************************
      subroutine inverse (dc,dci)
!**********************************************************************
!**********************************************************************
      double precision ::  dc(3,3),dci(3,3),det,deti
    
      det= dc(1,1)*( dc(2,2)*dc(3,3) - dc(3,2)*dc(2,3) ) &
        - dc(2,1)*( dc(1,2)*dc(3,3) - dc(3,2)*dc(1,3) )  &
        + dc(3,1)*( dc(1,2)*dc(2,3) - dc(2,2)*dc(1,3) )
      deti=1.d0/det

      dci(1,1)= (dc(2,2)*dc(3,3) - dc(3,2)*dc(2,3))*deti
      dci(1,2)=-(dc(1,2)*dc(3,3) - dc(3,2)*dc(1,3))*deti
      dci(1,3)= (dc(1,2)*dc(2,3) - dc(2,2)*dc(1,3))*deti

      dci(2,1)=-(dc(2,1)*dc(3,3) - dc(3,1)*dc(2,3))*deti
      dci(2,2)= (dc(1,1)*dc(3,3) - dc(3,1)*dc(1,3))*deti
      dci(2,3)=-(dc(1,1)*dc(2,3) - dc(2,1)*dc(1,3))*deti

      dci(3,1)= (dc(2,1)*dc(3,2) - dc(3,1)*dc(2,2))*deti
      dci(3,2)=-(dc(1,1)*dc(3,2) - dc(3,1)*dc(1,2))*deti
      dci(3,3)= (dc(1,1)*dc(2,2) - dc(2,1)*dc(1,2))*deti



      return
      end
!**********************************************************************
