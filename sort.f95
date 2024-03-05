! Code to make windows of randomly chosen dof value
! Made by Abhik Ghosh Moulick
! S N Bose National Centre for Basic Sciecnes
!Explanation of the parameters:
!maxframe=100000, maximum number of frames use for calculation set to 100000. Need to increase for larger trajectory.
!maxres=999,maximum number of entities (residues for protein, number of base pairs-1 for base pairs step parameters, number of base pairs for base pairs parameters, number of bases for backbone and sugar torsions in nucleic acid 
!maxwindow=20,maximum number of windows over which thermodynamics calculated
!maxconfpw=5000,maximum number of configuration per window, must be equal to maxframe/maxwindow
!jmax=2,considering the reference system and system of interest

integer, parameter :: maxframe=100000,maxres=999,maxwindow=20,maxconfpw=5000,jmax=2,iseed=187656
real*4, parameter:: tol=1.e-4 
character*1000 :: output_smp,input
character*50 :: char1,char2
character*50 :: resname(maxres)
character*50 apo_smp(maxwindow),apo_maxmin(maxwindow),apo_hist(maxwindow)
character*50 holo_smp(maxwindow),holo_maxmin(maxwindow),holo_hist(maxwindow)
!character*50 conftherm_win,conftherm_av
real*4 valuemean(maxwindow)
real*4 var(maxwindow)  
real*4 :: dofval(maxres,maxframe)
real*4 :: dfsample(maxwindow,maxconfpw)
!real*4 :: dof(maxframe,maxres,maxdf)
real*4 :: dof(maxframe,maxres)
!integer :: jmax

open(10,file='input_parameters.dat')

read(10,*) nframe !"Number of conformations you want to perform the analysis:the number must be below 100000"
read(10,*) nres !"Number of residues/nucleotides you want to perform the analysis:the number must be below 1000"
read(10,*)nwindow 
!read(10,*)nval !2*nwindow
read(10,*)nconf !"number of frames per window"
read(10,*)nbin ! 'number of bins in histogram"

call ranint(iseed)

!open(10,file='list')

do l=1,jmax

 
   read(10,*)input,output_smp   

 
   open(20,file=input)
   open(30,file='resname.dat')
	
   do i = 1,nframe
      do j = 1, nres		
	 read(20,*)nstep,resname(j),dof(i,j)              
      end do
   end do
   do j=1,nres
      write(30,*)j,resname(j)
   enddo
   close(20)
   close(30)
   

   do j = 1,nres        
      do i = 1,nframe
         dofval(j,i)=dof(i,j)
      end do
   end do
        

      do i=1,nres
         do k=1,nwindow
            valuemean(k) =0.
 	     open(40,file=output_smp)
            do j=1,nconf
               nsample= int(nframe*ranf())
               if (nsample.eq.0) nsample=1
               dfsample(k,j)=dofval(i,nsample)
               valuemean(k) = valuemean(k)+dfsample(k,j)
               write(40,*)j,resname(i),dfsample(k,j)
            enddo
            valuemean(k) = valuemean(k)/float(nconf)
            var(k) = 0.
            do j=1,nconf
               var(k) = var(k)+(dfsample(k,j)-valuemean(k))**2
            enddo
            var(k) = var(k)/float(nconf)
            var(k)= sqrt(var(k))   
         enddo

         iflag=0
         do k = 1,nwindow-1
           do kprime=k+1,nwindow
              diffmean = abs(valuemean(k)-valuemean(kprime))
              sumstd=.5*(var(k)+var(kprime))
              if(abs(diffmean-sumstd).gt.tol)iflag=1
           enddo
        enddo
     
        if(iflag.eq.0) then 
          write(*,*)'notconverged','syst',l,'res',resname(i)
          write(*,*)'PROGRAM will be terminated:needs longer trajectory'
          stop
        endif
    enddo          
  enddo


do i=1,nwindow
  read(10,*)apo_smp(i),apo_maxmin(i),apo_hist(i)
  read(10,*)holo_smp(i),holo_maxmin(i),holo_hist(i)
enddo
open(20,file='listhistprotmaxmin')
do i=1,nwindow
  write(20,*)apo_smp(i),apo_maxmin(i)
  write(20,*)holo_smp(i),holo_maxmin(i)
enddo
close(20)

open(20,file='listhistprot')
do i=1,nwindow
  write(20,*)apo_smp(i),apo_hist(i)
  write(20,*)holo_smp(i),holo_hist(i)
enddo
close(20)

open(20,file='histparameter')
write(20,*)2*nwindow
write(20,*)nconf
write(20,*)nbin
close(20)

open(20,file='parameter_maxmin')
write(20,*)2*nwindow
write(20,*)nconf
close(20)

open(20,file='num_window')
write(20,*)nwindow
close(20)

nbindash=nbin+2
open(20,file='listconformational_td')
write(20,*)nwindow,nbindash
!write(20,*)conftherm_win,conftherm_av
do i=1,nwindow
 write(20,*)holo_hist(i),apo_hist(i)
enddo
close(20)

close(10)
stop
end


      REAL function ranf()
!      common /rjran/ i3,i2,i1,i0

!     berkeley random number generator
!     range changed to 0 < 1

      INTEGER  I0, I1, I2, I3, J0, J1, J2, J3, K0, K1, K2, K3
      INTEGER  M0, M1, M2, M3, MM
      REAL     T1, T2, T3, T4     
      common /rjran/ i3,i2,i1,i0
      parameter (m3=647,m2=1442,m1=3707,m0= 373)
      parameter (t4=2.0**48,t3=2.0**36,t2=2.0**24,t1=2.0**12)
      parameter (mm=4096)


      ranf = float(i3)/t1 + float(i2)/t2 + float(i1)/t3 + float(i0)/t4
      if(ranf.ge.0.9999999)  ranf = 0.0

      j0 = m0 * i0
      j1 = m0 * i1 + m1 * i0
      j2 = m0 * i2 + m1 * i1 + m2 * i0
      j3 = m0 * i3 + m1 * i2 + m2 * i1 + m3 * i0
      k0 = j0
      k1 = j1 + k0 / mm
      k2 = j2 + k1 / mm
      k3 = j3 + k2 / mm
      i0 = mod(k0,mm)
      i1 = mod(k1,mm)
      i2 = mod(k2,mm)
      i3 = mod(k3,mm)
      return
      end      


      subroutine ranint(istart)

!    initialize random number generator

      INTEGER  IRAN(4), ISTART

      iran(1)=12345 + istart*1000
      iran(2)=12345 + istart*1000
      iran(3)=12345 + istart*1000
      iran(4)=12345 + istart*1000
      call ranset(iran)
      call ranget(iran)

      return
      end


      SUBROUTINE RANSET(IRAN)
 !     COMMON /RJRAN/ II3,II2,II1,II0

      INTEGER    IRAN(4), MM, NN
      PARAMETER (MM = 4096, NN = 100000)
      INTEGER    II0, II1, II2, II3
      INTEGER    I0, I1, I2, I3, J0, J1, J2, J3
      COMMON /RJRAN/ II3,II2,II1,II0
      i3 = iran(1)
      i2 = iran(2)
      i1 = iran(3)
      i0 = iran(4)
      call divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,mm,ii0)
      call divide(j3,j2,j1,j0,i3,i2,i1,i0,nn,mm,ii1)
      call divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,mm,ii2)
      call divide(j3,j2,j1,j0,i3,i2,i1,i0,nn,mm,ii3)
      return
      end


      subroutine ranget(iran)
!     common /rjran/ ii3,ii2,ii1,ii0

      INTEGER    mm, m10, m21, m20, m32, m31, m30
      common /rjran/ ii3,ii2,ii1,ii0
      parameter (mm = 100000)
      parameter (m10 = 4096)
      parameter (m21 =  167, m20 = 77216)
      parameter (m32 =    6, m31 = 87194, m30 = 76736)

      INTEGER    iran(4)
      INTEGER    ii3, ii2, ii1, ii0  
      INTEGER    J0, J1, J2, J3, K0, K1, K2, K3

      j0 = ii0 + m10 * ii1 + m20 * ii2 + m30 * ii3
      j1 =                   m21 * ii2 + m31 * ii3
      j2 =                               m32 * ii3
      j3 =                                       0
      k0 = j0
      k1 = j1 + k0 / mm
      k2 = j2 + k1 / mm
      k3 = j3 + k2 / mm
      iran(4) = mod(k0,mm)
      iran(3) = mod(k1,mm)
      iran(2) = mod(k2,mm)
      iran(1) = mod(k3,mm)
      return
      end


      subroutine divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,id,ir)

      INTEGER     I3,I2,I1,I0,J3,J2,J1,J0,ID,IR,NN,K0,K1,K2,K3

!     given the integer i = i0 + nn * (i1 + nn * (i2 + nn * i3))
!     this routine calculates j = i / id and ir = mod(i, id)
!     j is expressed as i, ir is just an integer

      j3 = i3 / id
      k3 = mod(i3, id)
      k2 = k3 * nn + i2
      j2 = k2 / id
      k2 = mod(k2, id)
      k1 = k2 * nn + i1
      j1 = k1 / id
      k1 = mod(k1, id)
      k0 = k1 * nn + i0
      j0 = k0 / id
      ir = mod(k0, id)
      return
      end
