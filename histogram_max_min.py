import numpy as np
import matplotlib.pyplot as plt
q1=[]
q2=[]
f = open('listhistprotmaxmin','r')
val = int(input("Enter your Number of data file in 'listhistprotmaxmin': "))
nconf= int(input("write Number of conformations in each file:"))
for j in range(val):
	s = f.readline()
	token = s.split()
	q1.append(token[0])
	q2.append(token[1])
	print(token[0])
	f1 = open(q1[j],'r')
	f2 = open(q2[j],'w')
	dof = np.zeros(nconf)
	for k in range(nconf):
		s = f1.readline()
		token = s.split()
		dof[k] = token[2]
		
	maximum = max(dof)
	minimum = min(dof)
	print(maximum,minimum,file=f2)
	
	
	
	
	
	
	
	
	
	
	
	


