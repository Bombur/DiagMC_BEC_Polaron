#!/usr/bin/env python3
import numpy as np
import math

pi=math.pi
relm=0.263158
qcs = [10,50,100,200,300,500,700, 1000,2000,3000,5000, 7000, 10000]
#alphas = [0.00001,0.0001,0.001,0.002,0.003, 0.005, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9, 1, 2,3,4,5]
alphas =[0.1,5]

def Eprenorm(alpha, qc):
	return alpha/np.sqrt(2)/pi *(1+ 1/relm)*qc

for a in qcs:
	for q  in alphas:
		print('%i \t %.g \t %.4g' %(a, q, Eprenorm(q, a)))
	print('#')
