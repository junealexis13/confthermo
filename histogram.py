import numpy as np
import matplotlib.pyplot as plt

f = open('listhistprot','r')
q1=[]
q2=[]
val = int(input("Enter your Number of data file in 'listhistprot': "))
nconf= int(input("write Number of conformations in each file:"))
nbin = int(input("write nbin:"))
for j in range(val):
	s = f.readline()
	token = s.split()
	q1.append(token[0])
	q2.append(token[1])
	print(token[0])
	f1 = open(q1[j],'r')
	f2 = open(q2[j],'w')
	dof = np.zeros(nconf)
	maximum,minimum = float(input("maximum:")),float(input("minimum:"))
	for k in range(nconf):
		s = f1.readline()
		token = s.split()
		dof[k] = token[2]
	binwidth = abs((maximum-minimum)/(nbin))
	print(binwidth, file=f2)
	xlist = np.linspace(minimum,maximum,nbin+1)
	hisdof = np.zeros(nbin+2)
	
	for k in range(nconf):
		a2 = dof[k] - minimum
		if (a2 % binwidth == 0.0):
			l2 = int(a2/binwidth)
		else: 
			l2 = int(a2/binwidth) + int(1)
		hisdof[l2] = hisdof[l2] + int(1)
	hisdof = hisdof/int(nconf)
	for m in range(nbin+2):
		print(m*binwidth+minimum,hisdof[m],file=f2)
	
	
	
	
	
	
	
	
	
	
	
	


