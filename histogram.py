import numpy as np
import matplotlib.pyplot as plt

f = open('listhistprot','r')
q1=[]
q2=[]
val = int(input("Enter your Number of data file in 'listhistprot': "))
npoint=int(input("write Number of conformations in each file:"))
binwidth = int(input("write binwidth:"))

for j in range(val):
	s = f.readline()
	token = s.split()
	q1.append(token[0])
	q2.append(token[1])
	print(token[0])
	f1 = open(q1[j],'r')
	f2 = open(q2[j],'w')
	dof = np.zeros(npoint)
	for k in range(npoint):
		s = f1.readline()
		token = s.split()
		dof[k] = token[2]
		
	maximum = max(dof)
	minimum = min(dof)
	binnumber = int(abs((maximum-minimum)/int(binwidth)))  # Binwidth
	print('Maximum value=',maximum,'minimum value=',minimum,"binwidth=",binwidth,"binnumber=",binnumber)	
	xlist = np.linspace(minimum,maximum,binnumber)
	hisdof = np.zeros(binnumber+2)
	
	for k in range(npoint):
		a2 = dof[k] - minimum
		if (a2 % binwidth == 0.0):
			l2 = int(a2/binwidth)
		else: 
			l2 = int(a2/binwidth) + int(1)
		hisdof[l2] = hisdof[l2] + int(1)
	hisdof = hisdof/int(npoint)
	for m in range(binnumber+2):
		print(m*binwidth+minimum,hisdof[m],file=f2)
	
	
	
	
	
	
	
	
	
	
	
	


