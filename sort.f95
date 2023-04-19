! Code to sort dihedral angle value for each residue.
! Made by Abhik Ghosh Moulick
! S N Bose National Centre for Basic Sciecnes
! This code is written to sort dihedral angle phi, psi, chi1 for protein
! If there are more than 3 conformational variable, one can need to add array for that.


integer,parameter :: frame=1000, nres=123  ! Change = Give frame number & residue number
character*1000 :: input1,input2,junk,output1
character*50 :: char1,char2,resname(500)
real*4 :: phi(frame,nres)
real*4 :: psi(frame,nres)
real*4 :: chi1(frame,nres)
real*4 :: chi2(frame,nres)
real*4 :: chi3(frame,nres)
real*4 :: chi4(frame,nres)
real*4 :: chi5(frame,nres)
integer :: jmax

write(*,*) "For how many system you want to perform the analysis, for two different system give 2, modify file 'list'  accordingly"
read(*,*) jmax

open(10,file='list')
do l=1,jmax
	do i = 1, frame
		do j = 1,nres
			phi(i,j) = 0.0
			psi(i,j) = 0.0
			chi1(i,j) = 0.0
			chi2(i,j) = 0.0
			chi3(i,j) = 0.0
			chi4(i,j) = 0.0
			chi5(i,j) = 0.0
		end do
	end do
	read(10,*)input1,output1
	write(*,*)input1
	open(20,file=input1)
	open(40,file=output1)
	do i = 1,frame
		do j = 1, nres
			read(20,*)nstep,resname(j),phi(i,j),psi(i,j),chi1(i,j),chi2(i,j),chi3(i,j),chi4(i,j),chi5(i,j) 
		end do
	end do
	
	do j = 1,nres
		do i = 1,frame
			write(40,4)j,resname(j),phi(i,j),psi(i,j),chi1(i,j),chi2(i,j),chi3(i,j),chi4(i,j),chi5(i,j)
4 			FORMAT (i4,1x,a4,1x,7(f10.4,1x))
		end do
	end do
end do
end
        
        
        
        
        
      
