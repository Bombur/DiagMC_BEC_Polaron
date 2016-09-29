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
peter_loc='/home/h/H.Guertner/theorie/BEC_Polaron/peter/'

#helper function
def find_nearest(array,value):
    return (np.abs(array-value)).argmin()
    

	

count=0
for folder in os.listdir(os.getcwd()):
	if os.path.exists(folder + "/fig4data"):
		count +=1
		
alpha = np.loadtxt('qc=200/fig4data')
alpha = alpha[:,0]
qcs = np.empty(count)

EpMC = np.empty([count, alpha.size])
EpMCerr = np.empty([count, alpha.size])
Ep = np.empty([count, alpha.size])
Eperr = np.empty([count, alpha.size])

it=0
for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/fig4data"):
		os.chdir(folder)
		#Import Parameters of each run
		with open('DiagMC_BEC.json') as parit_file:
			parit = json.load(parit_file)
			
		qcs[it] = parit["Q_Cutoff"]
		Eptmp = np.loadtxt('fig4data') 
		EpMC[it,:] = Eptmp[:,1]
		EpMCerr[it,:] = Eptmp[:,2]
		Eptmp = np.loadtxt('fig6data') 
		Ep[it,:] = Eptmp[:,1]
		Eperr[it,:] = Eptmp[:,2]
		
		os.chdir(os.pardir)
		it+=1

sort = np.argsort(qcs)
qcs = qcs[sort]
EpMC = EpMC[sort]
EpMCerr = EpMCerr[sort]
Ep = Ep[sort]
Eperr = Eperr[sort]


#fig5 plot
for plotit in range(alpha.size):
	plt.errorbar(qcs, Ep[:,plotit], Eperr[:, plotit], fmt='o', label='HP')

	plt.title('Alpha = %.1f' %(alpha[plotit]))
	plt.yscale('linear')
	plt.xlabel(r'$Q_c$')
	plt.ylabel(r'$Ep^{MC}+ E_{Ren}$')
	plt.xlim([0,6000])
	plt.savefig('Alpha%.1f' %(alpha[plotit]) + '.pdf')
	plt.clf()	
	
#fig4all
for plotit in range(count):
	plt.errorbar(alpha, EpMC[plotit,:], EpMCerr[plotit,:], label=r'$Q_c$ = ' + str(qcs[plotit]))

plt.yscale('log')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$-Ep^{MC}$')
plt.legend(loc=4)
plt.savefig('fig4all.pdf')
plt.clf()	

#fig6all
for plotit in range(count):
	plt.errorbar(alpha, Ep[plotit,:], Eperr[plotit,:], label=r'$Q_c$ = ' + str(qcs[plotit]))

plt.yscale('linear')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$Ep^{MC} + E_{Ren}$')
plt.legend(loc=3)
plt.savefig('fig6all.pdf')
plt.clf()	

#fig6small
for plotit in range(count):
	plt.errorbar(alpha, Ep[plotit,:], Eperr[plotit,:], label=r'$Q_c$ = ' + str(qcs[plotit]))

Gr = np.genfromtxt(peter_loc + 'Grsmall_Vlie.csv',delimiter=',', skip_header=6)
Grerr = np.genfromtxt(peter_loc +'Grsmall_Vlie_err.csv',delimiter=',', skip_header=6)
Grerr[:,1]-=Gr[:,1]
FeyVar = np.genfromtxt(peter_loc + 'Grsmall_FM.csv',delimiter=',', skip_header=6)
MF = np.genfromtxt(peter_loc + 'Grsmall_MF.csv',delimiter=',', skip_header=6)

plt.errorbar(Gr[:, 0], Gr[:, 1], Grerr[:,1], label="Vlietinck in Grusdt paper")
plt.plot(FeyVar[:, 0], FeyVar[:, 1], label = "FeyVar")
plt.plot(MF[:, 0], MF[:, 1], label = "MF")

plt.yscale('linear')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$Ep^{MC} + E_{Ren}$')
#plt.legend(loc=3)
plt.savefig('fig6small.pdf')
plt.clf()	

