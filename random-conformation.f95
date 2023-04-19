! Code to get random conformation from dihedral angle dataset.
! Made by Abhik Ghosh Moulick
! S N Bose National Centre for Basic Sciecnes
! This code is written to choose random residue and conformations from dihedral angle dataset
! This data can be used for intial convergence of the system 
! If there are more than 3 conformational variable, one could add array for that.


integer,parameter :: frame=1000, nres=123  ! Change = change accordingly frame number and number of residue
character*1000 :: input1,input2,junk,output1,output2
character*50 :: char1,char2,resname(nres)  
real*4 :: phi(nres,frame)
real*4 :: psi(nres,frame)
real*4 :: chi1(nres,frame)
real*4 :: chi2(nres,frame)
real*4 :: chi3(nres,frame)
real*4 :: chi4(nres,frame)
real*4 :: chi5(nres,frame)
integer :: seed, num_res, num_frame, k, p, num_dof
real :: a
integer :: jmax

write(*,*) "For how many system you want to perform the analysis, for two different system give 2, modify file 'list'  accordingly"
read(*,*) jmax

write(*,*) "How many conformation you want to randomly generate"
read(*,*) p 

seed = 187656
open(10,file='list')
do l=1,jmax
	do i = 1, nres
		do j = 1,frame
			phi(i,j) = 0.0
			psi(i,j) = 0.0
			chi1(i,j) = 0.0
			chi2(i,j) = 0.0
			chi3(i,j) = 0.0
			chi4(i,j) = 0.0
			chi5(i,j) = 0.0
		end do
	end do
        read(10,*)input1,output1,output2
        open(20,file=output1)
        open(40,file=output2)
        do j = 1,nres
        	do i = 1, frame
        		read(20,*)nstep,resname(j),phi(j,i),psi(j,i),chi1(j,i),chi2(j,i),chi3(j,i),chi4(j,i),chi5(j,i)
		end do
	end do
	
!       Random generation of conformation
	
	call system_clock(k)
  	call srand(seed)
  	write(*,*) "For which residue you want to perform conformational thermodynamics, give residue number"
  	read(*,*) num_res
  	write(*,*)"Residue of interest is ",resname(num_res)
  	write(*,*)"For which dihedfral you want to perform conformational thermodynamics"
  	write(*,*)"write 1 for phi"
  	write(*,*)"write 2 for psi"
  	write(*,*)"write 3 for chi1"
  	write(*,*)"write 4 for chi2"
  	write(*,*)"write 5 for chi3"
  	write(*,*)"write 6 for chi4"
  	write(*,*)"write 7 for chi5" 
  	read(*,*)num_dof
	do i = 1, p
  		a = 1000*rand(seed)
  		seed = seed + k
  		num_frame = int(a)
  		if (num_dof .eq. 1) then
			write(40,50)i,resname(num_res),phi(num_res,num_frame)
		else if (num_dof .eq. 2) then
			write(40,50)i,resname(num_res),psi(num_res,num_frame)
		else if (num_dof .eq. 3) then
			write(40,50)i,resname(num_res),chi1(num_res,num_frame)
		else if (num_dof .eq. 4) then
			write(40,50)i,resname(num_res),chi2(num_res,num_frame)
		else if (num_dof .eq. 5) then
			write(40,50)i,resname(num_res),chi3(num_res,num_frame)
		else if (num_dof .eq. 6) then
			write(40,50)i,resname(num_res),chi4(num_res,num_frame)
		else if (num_dof .eq. 7) then
			write(40,50)i,resname(num_res),chi5(num_res,num_frame)
		end if
50		format(i4,1x,a4,1x,1(f10.4,1x))
 	end do
end do
end
        
        
        
        
        
      
