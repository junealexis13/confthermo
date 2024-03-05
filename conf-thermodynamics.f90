
        parameter (maxbin=100,maxwindow=100)           !change
        common histdof(maxbin),histdof_c(maxbin),histdof_f(maxbin)
        real histmaxdof_c(maxbin)
        integer maxdof_c(maxbin),maxdof_f(maxbin)
        real histmaxdof_f(maxbin)
        real curvdof_f(maxbin)
 
        real curvdof_c(maxbin)

        dimension free_enDiff(maxwindow),delS(maxwindow)
       
        character*1000 input1,input2
         
 
       open(90,file='listconformational_td')
       read(90,*)nwindow,nbin
 
       
       open(50,file='conftherm_window.dat')
       open(60,file='conftherm_av.dat')
       
       do l=1,nwindow
            
          read(90,*)input1,input2
            
             
	
!*****************************************************************
  	  
!   for the complex

!******** Read the input histogram file ******************
	  
	  
!	  
	  open(10,file=input1)
	  	  
          read(10,*)width  
	  
	  
    	  do j=1,nbin
  	     read(10,*)ang,histdof(j)
 		
  	     if(histdof(j).eq.0.0)histdof(j)=1.0e-06
  	     histdof_c(j)=histdof(j)
 ! 					
	     enddo
	     	      
	     entropy_c=entropy(nbin,width)
	     
	call histcal(histdof_c,curvdof_c,histmaxdof_c,maxdof_c,maxndof_c,nbin,width)
!		
!*****************************************************************
  	  
!   for the free protein

!***************************************************************** 
!******** Read the input histogram file ******************
     	    open(20,file=input2)

           read(20,*)width

    	   do j=1,nbin
    	      
	      read(20,*)ang,histdof(j)  
	      if(histdof(j).eq.0.0)histdof(j)=1.0e-06
	      histdof_f(j)=histdof(j)
	         
	   enddo

           entropy_f=entropy(nbin,width)
	call histcal(histdof_f,curvdof_f,histmaxdof_f,maxdof_f,maxndof_f,nbin,width)

        Fdof= Fdihed(histmaxdof_f,curvdof_f,histmaxdof_c,curvdof_c,maxndof_f,maxndof_c,width)

!*****************************************************************
    	   delS(l) = entropy_c-entropy_f
    	    
          free_enDiff(l)= Fdof
           
          close(10)
          close(20)
           
          write(50,*)l,2.5*free_enDiff(l),2.5*delS(l)

        enddo
   
    	
    	aven=0.
        do l = 1,nwindow
    	   aven = aven+free_enDiff(l)
    	enddo
    	aven = aven/float(nwindow)
    	error=0.
    	do l = 1,nwindow
    	   error=error+(free_enDiff(l)-aven)**2
    	enddo
    	error = sqrt(error/float(nwindow))
    	erroren = .5*error/sqrt(float(nwindow))
    	
    	avent=0.
        do l = 1,nwindow
    	   avent = avent+delS(l)
    	enddo
    	avent = avent/float(nwindow)
    	error=0.
    	do l = 1,nwindow
    	   error=error+(delS(l)-avent)**2
    	enddo
    	error = sqrt(error/float(nwindow))
    	errorent = .5*error/sqrt(float(nwindow))
    	
    	write(60,*)'aven=',aven,'erroren=',erroren,'avent=',avent,'errorent=',errorent
    	close(50)
    	close(60)
    	  
        stop
	end



        function entropy(nbin,width) 
        parameter (maxbin=100)           !change
        common histdof(maxbin),histdof_c(maxbin),histdof_f(maxbin)
        
        
	
         value1= histdof(1)*log(histdof(1))
            
         valueN= histdof(nbin)*log(histdof(nbin))
            

           sum = 0.5*(value1+valueN)
          
            do k=2,nbin-1
               
               xint=histdof(k)*log(histdof(k))
                
               sum = sum + xint
            
            enddo
            
            entropy=-sum*width
            
            return
            end


!************************************************************************************ 
  	 !************************************************************************************ 
  	  
!   routine for calculation of maximum of a histogram and curvature near the maximum  	  

!************************************************************************************ 
	
	subroutine histcal(hist,curv,histmax,maxv,maxn,nbin,h)
	parameter(maxbin=100)

	integer s,maxn,bigl(maxbin)
        real hist(maxbin)
        real curv(maxbin),histmax(maxbin),big(maxbin),x1,x2
        real dum(maxbin),h
	integer j,k,maxv(maxbin)
	

  	  do j=1,nbin
	    dum(j) = hist(j)	  
	  enddo
	  call maxcal(dum,nbin,big,bigl,s)
	 maxn=s
	  
	  do k=1,s
	    histmax(k) = big(k)
	    maxv(k) = bigl(k)
	  enddo

!*****************************************************************
!	to calculate the curvature near maximum of histogram
!*****************************************************************

          do k=1,s
            if(histmax(k) .ne. 0.0) then
              if(dum(maxv(k)+1)-dum(maxv(k)) .ne. 0.0 .and. dum(maxv(k)-1)-dum(maxv(k)) .ne. 0.0) then
                x1 = (-1.0)*2.0*(dum(maxv(k)+1)-dum(maxv(k)))/(h*h)
                x2 = (-1.0)*2.0*(dum(maxv(k)-1)-dum(maxv(k)))/(h*h)
                curv(k) = (x1+x2)/2.0

              else if(dum(maxv(k)+1)-dum(maxv(k)) .eq. 0.0) then
                curv(k) = (-1.0)*2.0*(dum(maxv(k)-1)-dum(maxv(k)))/(h*h)              

              else if(dum(maxv(k)-1)-dum(maxv(k)) .eq. 0.0) then
                curv(k) = (-1.0)*2.0*(dum(maxv(k)+1)-dum(maxv(k)))/(h*h)

              else 

              endif
            else 
              curv(k) = 0.0
            endif

          enddo
   
		
          end

!************************************************************************************ 
  	  
!   routine for calculation of local maxima of an array  	  

!************************************************************************************ 

	subroutine maxcal(x,n,maxv,maxl,m)
      parameter(maxbin=100)
	integer i,j,n,m,maxl(maxbin),sloc
	real x(n),maxv(maxbin),maxm

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
	
	  
1	m=j

	end
  
 
	!************************************************************************************ 
  	  
!   routine for calculation of free energy and entropy contribution of any dna dihedral  	  

!************************************************************************************ 

	
	Function Fdihed(histmaxd_f,curvd_f,histmaxd_c,curvd_c,maxn_f,maxn_c,h)
!	implicit none
	parameter(maxbin=100,pi=3.141)
	integer nres,i,j,k,maxn_f, maxn_c,maxv_f,maxv_c
	real histmaxd_f(maxbin),curvd_f(maxbin),histmaxd_c(maxbin),curvd_c(maxbin)
	real Fdihed,histmax_f(maxbin),histmax_c(maxbin)
	real F,S,w
	real hsum_f,hsum_c,hcorrf,hcorrc,h,Fcorr,Fmin

	
	
	
	  hsum_f = 0.0
	  
	  do j=1,maxn_f
	    hsum_f = hsum_f + histmaxd_f(j)
	  enddo

	  hsum_c = 0.0
	  do j=1,maxn_c
	    hsum_c = hsum_c + histmaxd_c(j)
	  enddo	  
	  
	  Fdihed = 0.0
	   
	  do j=1,maxn_f
	    do k=1,maxn_c
	      if(maxn_f.eq.1.and.maxn_c.eq.1) then 
	        w = 1
	      else if(maxn_f.eq.2.and.maxn_c.eq.1) then 
	        w = histmaxd_f(j)/hsum_f
	      else if(maxn_f.eq.1.and.maxn_c.eq.2) then 
	        w = histmaxd_c(k)/hsum_c
	      else 
	        w = histmaxd_f(j)*histmaxd_c(k)/(hsum_f + hsum_c)
	      endif
	      
	      hcorrf = sqrt(2.0*pi*histmaxd_f(j)/curvd_f(j))
	      hcorrc = sqrt(2.0*pi*histmaxd_c(k)/curvd_c(k))

	      Fcorr = -(hcorrc - hcorrf)
	      F = -log(histmaxd_c(k)/histmaxd_f(j)) + Fcorr
	      Fmin = -log(histmaxd_c(k)/histmaxd_f(j))

	      
	      Fdihed = Fdihed + w*F

	    enddo
	  enddo	      
	
	end

