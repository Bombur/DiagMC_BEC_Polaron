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
    
readin = ['fig4_qc100_MC.csv', 'fig4_qc200_MC.csv', 'fig4_qc3000_MC.csv']
folders = ['qc=100','qc=200','qc=3000']

alpha = [1,2,3,4,5]
alpha = np.asarray(alpha)
qcs = np.empty(3)

EpMC = np.empty([qcs.size,alpha.size] )
EpMCerr = np.empty([qcs.size, alpha.size])
Vlie = np.empty([qcs.size, alpha.size])

it=0
for folder in folders:
	os.chdir(folder)
	with open('DiagMC_BEC.json') as parit_file:
		parit = json.load(parit_file)
		
	qcs[it] = parit["Q_Cutoff"]
	Eptmp = np.loadtxt('fig4data') 
	EpMC[it,:] = Eptmp[4:,1]
	EpMCerr[it,:] = Eptmp[4:,2]
	Eptmp = np.genfromtxt(parit["Peter_Path"] + '/' + readin[it], delimiter=',', skip_header=6) 
	for line in Eptmp:
		if line[0] > 5.5:
			continue
		idx = find_nearest(alpha, line[0])
		Vlie[it,idx] = line[1]
			
	os.chdir(os.pardir)
	it+=1

for qit in range(qcs.size):
	plt.errorbar(alpha, (EpMC[qit, :]/Vlie[qit,:]), (EpMCerr[qit,:]/Vlie[qit,:]), fmt= 'o', label=r'$Q_c = %.i' %(qcs[qit]))
	plt.axhline(y = np.sqrt(2), color ='g')

plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\frac{E_p^{MC}}{E_p^{Vlietinck}}$')
plt.legend(loc=4)
plt.savefig('fig4_comp_Vlie.pdf')
plt.clf()

