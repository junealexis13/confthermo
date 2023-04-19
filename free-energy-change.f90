		implicit none
		integer nresprot,nresdna,maxbin,nn,numbin
        parameter (nresprot=1,numbin=10000,nn=100)					!!!!!!!!!!!change
        real histdof(nresprot,numbin),pi
        real histmaxdof_f(500,20),delF(nresprot)
        real curvdof_f(500,20)
        real histmaxdof_c(500,20), delFtotp
        real curvdof_c(500,20)
        integer maxndof_f(500)
        integer maxndof_c(500)
        real delFbind1,delSbind1,delFbind2,delSbind2,delFbind3,delSbind3       

        real angle(numbin),h
		integer i,j,s,l,m,junk,resnump(500),resnumd(500),res(500)	
        integer maxdof_f(500,20)
        integer maxdof_c(500,20), nline1
		character*20 resdp(500),resdd(500),junk2
		real Fdof(nresprot)
		real Sdof(nresprot),delStotp,delS(nresprot)
        real delFdof
        real delSdof
		character*1000 input1,input2,output
        real delFbr(100),delSbr(100),delFbrtot,delSbrtot
        integer nbr,jbr,resbr1(100),resbr2(100),jmax

	data pi/3.14/	
	write(*,*)"write Number of file given in 'list-conformational-td' file"
    read(*,*)jmax
    
    write(*,*)"write the binwidth of your dataset"
    read(*,*)h
	h = h*(pi/180.0) 
	write(*,*) h
	open(90,file='list-conformational-td')
	do l=1,jmax
            read(90,*)input1,input2,output
            write(*,*)input1
	
!*****************************************************************
  	  
!   for the complex

!***************************************************************** 

!******** Count number of bin in the file ***************

        open(10,file=input1)	
        nline1 = 0
		DO
			READ (10,*, END=100)
			nline1 = nline1 + 1
		END DO
100 		CLOSE (10)
        
        maxbin = nline1 
        !write(*,*) maxbin
!******** Read the input histogram file ******************
		open(10,file=input1)	
		do i=1,nresprot
    	  do j=1,maxbin
	    read(10,*)angle(j),histdof(i,j)
	    !write(*,*)angle(j),histdof(i,j)
		enddo	  
		enddo
	close(10)
                	
	call histcal(histdof,curvdof_c,histmaxdof_c,maxdof_c,nresprot,maxndof_c,maxbin,h)
	
			
		
!*****************************************************************
  	  
!   for the free protein

!***************************************************************** 

!******** Count number of bin in the file ***************
		write(*,*)input2
        open(20,file=input2)	
        nline1 = 0
		DO
			READ (20,*, END=99)
			nline1 = nline1 + 1
		END DO
99 		CLOSE (20)
        
        maxbin = nline1 
        !write(*,*) maxbin
!******** Read the input histogram file ******************
	open(20,file=input2)
	do i=1,nresprot
    	  do j=1,maxbin
	    read(20,*)angle(j),histdof(i,j)
	  enddo
	enddo
	close(20)

	call histcal(histdof,curvdof_f,histmaxdof_f,maxdof_f,nresprot,maxndof_f,maxbin,h)
	
		
!*****************************************************************

!   calculation of Fs: protein

        call Fresdihed(histmaxdof_f,curvdof_f,histmaxdof_c,curvdof_c,nresprot,maxndof_f,maxndof_c,Fdof,Sdof,h)
        
	
!*****************************************************************



!	total free energy cost

        open(50,file=output)                 !change
  
    	delFtotp = 0.0
    	delStotp = 0.0  

        delFdof = 0.0
        delSdof = 0.0
        

        jbr=1

	do i=1,nresprot
	  delF(i) = Fdof(i)       !!total free energy for each residue
	  delS(i) = Sdof(i)      !!total entropy for each residue
   	  write(50,1)Fdof(i)*2.52
   	    
1	  FORMAT (2(f10.5,'  '))  

    enddo

    	
    end do  

    close(50)
    end
    	
    	  


!************************************************************************************ 
  	  
!   routine for calculation of maximum of a histogram and curvature near the maximum  	  

!************************************************************************************ 

	subroutine histcal(hist,curv,histmax,maxv,nres,maxn,maxbin,h)
	
	implicit none
	integer nres,maxbin,s,maxn(500),bigl(20)
        real hist(nres,maxbin)
        real curv(500,20),histmax(500,20),big(20),x1,x2
        real dum(maxbin),h,dum1(maxbin)
	integer i,j,k,maxv(500,20)
	
	do i=1,nres
  	  do j=1,maxbin
	    dum(j) = hist(i,j)	  
	  enddo
	  call maxcal(dum,maxbin,big,bigl,s)
	  maxn(i) = s
	  
	  do k=1,maxn(i)
	    histmax(i,k) = big(k)
	    maxv(i,k) = bigl(k)
	  enddo

!*****************************************************************
!	to calculate the curvature near maximum of histogram
!*****************************************************************

          do k=1,maxn(i)
            if(histmax(i,k) .ne. 0.0) then
              if(dum(maxv(i,k)+1)-dum(maxv(i,k)) .ne. 0.0 .and. dum(maxv(i,k)-1)-dum(maxv(i,k)) .ne. 0.0) then
                x1 = (-1.0)*2.0*(dum(maxv(i,k)+1)-dum(maxv(i,k)))/(h*h)
                x2 = (-1.0)*2.0*(dum(maxv(i,k)-1)-dum(maxv(i,k)))/(h*h)
                curv(i,k) = (x1+x2)/2.0
!                write(*,*)'here'
              else if(dum(maxv(i,k)+1)-dum(maxv(i,k)) .eq. 0.0) then
                curv(i,k) = (-1.0)*2.0*(dum(maxv(i,k)-1)-dum(maxv(i,k)))/(h*h)              
!                write(*,*)'hi'
              else if(dum(maxv(i,k)-1)-dum(maxv(i,k)) .eq. 0.0) then
                curv(i,k) = (-1.0)*2.0*(dum(maxv(i,k)+1)-dum(maxv(i,k)))/(h*h)
!                write(*,*)'bye'
              else 
!                write(*,*)'ho'
              endif
            else 
              curv(i,k) = 0.0
            endif
!            write(*,*)curv(i,k),histmax(i,k),dum(maxv(i,k)+1),dum(maxv(i,k)),maxv(i,k)
          enddo
        enddo
		
	end

!************************************************************************************ 
  	  
!   routine for calculation of local maxima of an array  	  

!************************************************************************************ 

	subroutine maxcal(x,n,maxv,maxl,m)

	integer i,j,n,m,maxl(20),sloc
	real x(n),maxv(20),maxm

	j=0
	do i=2,n-1
	  if(x(i).gt.x(i-1).and.x(i).gt.x(i+1).and.x(i).ge.0.0) then 
	    j=j+1
	    maxl(j)=i
	    maxv(j)=x(i)
	  else if(x(i).gt.x(i-1).and.x(i).ge.x(i+1).and.x(i).ge.0.0) then 
	   if(x(i).gt.x(i+2)) then
	    j=j+1
	    maxl(j)=i
	    maxv(j)=x(i)
	   endif
	  endif	  	  
	enddo
	
!	call maxov(x,n,maxm,sloc)
	
!	do i=1,j
!	  if(maxv(i).eq.maxm)goto 1
!	enddo
!	write(*,*)'here',j
!	j=j+1
!	maxl(j+1)=sloc
!	maxv(j+1)=maxm
	  
1	m=j

	end
  
!************************************************************************************ 
  	  
!   routine for calculation of global maximum of an array  	  

!************************************************************************************ 

	subroutine maxov(x,n,maxv,q)

	integer i,n,q
	real x(n),maxv,maxm

	maxm=x(1)
	do i=1,n
	  if(x(i).ge.maxm) then
	    maxm=x(i) 
	    q=i
	  endif
	enddo
	maxv=maxm

	end

!************************************************************************************ 
  	  
!   routine for calculation of free energy and entropy contribution of any dna dihedral  	  

!************************************************************************************ 

	subroutine Fresdihed(histmaxd_f,curvd_f,histmaxd_c,curvd_c,nres,maxn_f,maxn_c,Fdihed,Sdihed,h)
	
	implicit none
	integer nres,i,j,k,maxn_f(500),maxn_c(500),maxv_f,maxv_c,maxbin
	real histmaxd_f(500,20),curvd_f(500,20),histmaxd_c(500,20),curvd_c(500,20)
	real Fdihed(nres),Sdihed(nres),histmax_f(20),histmax_c(20)
	real F,S,w
	real hsum_f,hsum_c,hcorrf,hcorrc,pi,h,Fcorr,Fmin

	
	
	do i=1,nres
	  hsum_f = 0.0
	  do j=1,maxn_f(i)
	    hsum_f = hsum_f + histmaxd_f(i,j)
	  enddo
	  hsum_c = 0.0
	  do j=1,maxn_c(i)
	    hsum_c = hsum_c + histmaxd_c(i,j)
	  enddo	  
	  
	  Fdihed(i) = 0.0
	  Sdihed(i) = 0.0
	  do j=1,maxn_f(i)
	    do k=1,maxn_c(i)
	      if(maxn_f(i).eq.1.and.maxn_c(i).eq.1) then 
	        w = 1
	      else if(maxn_f(i).eq.2.and.maxn_c(i).eq.1) then 
	        w = histmaxd_f(i,j)/hsum_f
	      else if(maxn_f(i).eq.1.and.maxn_c(i).eq.2) then 
	        w = histmaxd_c(i,k)/hsum_c
	      else 
	        w = histmaxd_f(i,j)*histmaxd_c(i,k)/(hsum_f + hsum_c)
	      endif
	      hcorrf = sqrt(2.0*pi*histmaxd_f(i,j)/curvd_f(i,j))
	      hcorrc = sqrt(2.0*pi*histmaxd_c(i,k)/curvd_c(i,k))
!	      F = -log(histmaxd_c(i,k)/histmaxd_f(i,j))
	      Fcorr = -(hcorrc - hcorrf)
	      F = -log(histmaxd_c(i,k)/histmaxd_f(i,j)) + Fcorr
	      Fmin = -log(histmaxd_c(i,k)/histmaxd_f(i,j))
!	      write(*,7)Fmin,(-Fcorr),F,hcorrf,hcorrc,-curvd_c(i,k),-curvd_f(i,j)
7	      format(7(f8.4,'	'))
	      S = 0.5*log(curvd_f(i,j)/curvd_c(i,k))	      
	      Fdihed(i) = Fdihed(i) + w*F
	      Sdihed(i) = Sdihed(i) + w*S
!	      write(*,*)w,maxn_f(i),maxn_c(i),histmaxd_f(i,j),histmaxd_c(i,k),i
	    enddo
	  enddo	      
	enddo

	end



