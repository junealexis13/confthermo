# This code can be used to calculate mean and error
#import library

import matplotlib.pyplot as plt
import numpy as np
import math
import random
import statistics


f1 = open('all-freeenergychange.dat','r')
f2 = open('all-entropycost.dat','r')
f3 = open('final-value.dat','w')

datafile = int(input("Write the number of datapoint, for which average has been performed: "))

val_energy = np.zeros(datafile)
val_entropy = np.zeros(datafile)

for k in range(datafile):
	s = f1.readline()
	token = s.split()
	val_energy[k] = token[0]
	s = f2.readline()
	token = s.split()
	val_entropy[k] = token[0]

sd_energy = statistics.stdev(val_energy)
err_energy = sd_energy/(2*math.sqrt(datafile))

sd_entropy = statistics.stdev(val_entropy)
err_entropy = sd_entropy/(2*math.sqrt(datafile))
print('Final value with error is printed in the file "final-value.dat"')
print('average free energy',statistics.mean(val_energy),'	error in energy',err_energy,file=f3)
print('average entropy',statistics.mean(val_entropy),'	error in entropy',err_entropy,file=f3)


	
