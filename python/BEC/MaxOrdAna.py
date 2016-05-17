#!/usr/bin/env python3
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, errno, sys
import math
import json
#from scipy.optimize import curve_fit
from scipy import optimize
from scipy.interpolate import interp1d, splev, splrep


# definitions
pi = math.pi
scale= 'log'

#helper function
def find_nearest(array,value):
    return (np.abs(array-value)).argmin()
    	

count=0
for folder in os.listdir(os.getcwd()):
	if os.path.exists(folder + "/Epconvord"):
		count +=1
		alpha = np.loadtxt(folder + '/fig4data')
		
alpha = alpha[:,0]
qcs = np.empty(count)

ConvOrd = np.zeros((count, alpha.size))
MOOrd = np.zeros((count, alpha.size))

it=0
for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/fig4data"):
		os.chdir(folder)
		#Import Parameters of each run
		with open('DiagMC_BEC.json') as parit_file:
			parit = json.load(parit_file)
			
		qcs[it] = parit["Q_Cutoff"]
		MOtmp = np.loadtxt('Epconvord') 
		rows,cols = MOtmp.shape
		for i in range(rows):
			idx = find_nearest(alpha, MOtmp[i,0])
			ConvOrd[it, idx] = MOtmp[i,2] if MOtmp[i,2] > ConvOrd[it, idx] else ConvOrd[it, idx]
			MOOrd[it, idx] = MOtmp[i,3] if MOtmp[i,3] > MOOrd[it, idx] else MOOrd[it, idx]
		os.chdir(os.pardir)
		it+=1

sort = np.argsort(qcs)
qcs = qcs[sort]
ConvOrd = ConvOrd[sort]
MOOrd = MOOrd[sort]


for plotit in range(alpha.size):
	plt.plot(qcs, ConvOrd[:,plotit], 'o', label='Alpha = %.1f' %(alpha[plotit]))

plt.legend()
plt.yscale('linear')
plt.xlabel(r'$Q_c$')
plt.ylabel('Maximum Order relevant for Estimator')
plt.xlim([0,6000])
plt.savefig('ConvOrd.pdf')
plt.clf()	
	
for plotit in range(alpha.size):
	plt.plot(qcs, MOOrd[:,plotit], 'o', label='Alpha = %.1f' %(alpha[plotit]))

plt.legend()
plt.yscale('linear')
plt.xlabel(r'$Q_c$')
plt.ylabel('Maximum Order reached')
plt.xlim([0,6000])
plt.savefig('MOges.pdf')
plt.clf()	
