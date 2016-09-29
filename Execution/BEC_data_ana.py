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
from multiprocessing import Pool
sys.path.append("/project/theorie/h/H.Guertner/BEC_Polaron/python")
import Trafos
import Ep_ana as Epa

# definitions
pi = math.pi
scale= 'linear'
SEscale ='log'

# import parameters
with open('DiagMC_BEC.json') as params_file:
    params = json.load(params_file)

taus = params['Taus']
taumax = taus[-1]
bins = params['Bins']
peter_loc = params["Peter_Path"]
p  = params['Momentum']
#Loading Ws
data = np.loadtxt('data/Ep/Epvsws')
ws = data[:,0]
Eprenorm = params['Alpha']/np.sqrt(2)/pi *(1+ 1/params['Impurity_Mass'])*params['Q_Cutoff']
smend = 100 #plot limit for small plots

#helper function
def find_nearest(array,value):
    return (np.abs(array-value)).argmin()
    
        
Epa.mat_scripts()    #Mathematica Scripts

Epa.Stat_ana()		#Statistics Analysis


#---------------- plot Green-function

#sample
data = np.loadtxt('data/all_orders')
nozero = np.nonzero(data[:,1])
data = data[nozero]

#small plot
plt.errorbar(data[:smend,0],data[:smend,1], data[:smend,2], label="Sampled")
plt.yscale('log')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G(\tau)$')
plt.legend(loc = 4)
plt.savefig('ana/plot_all_small.pdf')
plt.clf()


#Big Plot and analysis
fitFile = open("ana/PolaronEnergy", "w")
fitFile.write("#Polaron Energy \n")
fitFile.write("Ep \t Ep Error \t Z \t Z Error \n")
fitFile.close()

#fitting	
nozero = np.nonzero(data[:,2])
fitdata = data[nozero]
	
fitstart = 100
fitend = -5
	
#Ep, Eperr, Z, _ = Epa.Ep_fit(fitstart, fitend, fitdata[:,0], fitdata[:,1], fitdata[:,2])
#if Ep !=0:
#	plt.plot(data[:,0], Epa.Gtaugrwp(data[:,0], Ep, Z),  label ='Ep = %.g +- %.e' %(Ep, Eperr))
	
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], fmt = 'bx', label='Sampled')		

plt.yscale('log')
plt.xlabel(r'$\tau$')
if params['Self_Energy']:
	plt.ylabel(r'$G_0 \Sigma(\tau)$')
else:
	plt.ylabel(r'$G(\tau)$')
plt.legend(loc = 3)
plt.savefig('ana/plot_all.pdf')
plt.clf()


#-------------- plot fake-function
#sample
fake = np.loadtxt('data/zero_order')
nozero = np.nonzero(fake[:,1])
fake = fake[nozero]
plt.errorbar(fake[:, 0], fake[:, 1], fake[:, 2], label='Sampled')

# analytical
def fake_integrand(tau):
  disp = np.power(params['Momentum'], 2)/2. - params['Chemical_Potential']
  return np.exp(- disp * tau)

fake_integrand = np.vectorize(fake_integrand)
plt.plot(fake[:, 0], fake_integrand(fake[:, 0]), label='Analytical result')
plt.yscale('linear')
plt.xlabel(r'$\tau$')
plt.ylabel(r'Fake-function')
plt.legend()
plt.savefig('ana/plot_0.pdf')
plt.clf()

#small plot
plt.errorbar(fake[:smend,0],fake[:smend,1], fake[:smend,2], label="Sampled")
plt.plot(fake[:smend, 0], fake_integrand(fake[:smend, 0]), label='Analytical result')
plt.yscale('linear')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G(\tau)$')
plt.legend()
plt.savefig('ana/plot_0_small.pdf')
plt.clf()


#--------------- plot 1st order
#sample
data = np.loadtxt('data/first_order')
datasmall = data
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')

#Mathematica

if params["FOG0SEG0"]:
	if not os.path.exists("./data/mat_1st_g0seg0"):
		os.system("math -script 1st_G0SEG0.m")
	try:
		mat=np.loadtxt("data/mat_1st_g0seg0")
		plt.plot(mat[:, 0], mat[:,1], label='Mathematica p=0')	
	except:
		pass
else:
	if not os.path.exists("./data/mat_1st_g0se"):
		os.system("math -script 1st_G0SE.m")
	try:
		mat=np.loadtxt("data/mat_1st_g0se")
		plt.plot(mat[:, 0], mat[:,1], label='Mathematica p=0')	
	except:
		pass

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_0 \Sigma G_0 (0, \tau)$')
plt.legend()
plt.savefig('ana/plot_1.pdf')
plt.clf()

#small plot
plt.errorbar(datasmall[:smend,0],datasmall[:smend,1], datasmall[:smend,2], label="Sampled")
try:
	plt.plot(mat[:smend, 0], mat[:smend, 1], label='Mathematica')
except:
	pass
plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G(\tau)$')
plt.legend()
plt.savefig('ana/plot_1_small.pdf')
plt.clf()



#--------------- plot 2nd order a
#sample
data = np.loadtxt('data/second_ordera')
datasmall= data
datacomp= data
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled 2a')
data2 = np.loadtxt('data/second_orderb')

"""
if params["Self_Energy"]:
	datacomp[:,1:] += data2[:,1:]
	datacomp= datacomp[nozero]
	plt.errorbar(datacomp[:, 0], datacomp[:, 1], datacomp[:,2], label='Complete 2')

	#if not os.path.exists("./data/mat_2nd_g0se_complete"):
	#	os.system("math -script 2nd_G0SE_Complete.m")
	#try:
	#	mat=np.loadtxt("data/mat_2nd_g0se_complete")
	#	plt.plot(mat[:, 0], mat[:,1], label='Mathematica p=0 Complete')	
	#except:
	#	pass
	
	#if not os.path.exists("./data/mat_2nd_g0se"):
	#	os.system("math -script 2nd_G0SE_Rainbow.m")
	#try:
	#	mat=np.loadtxt("data/mat_2nd_g0se")
	#	plt.plot(mat[:, 0], mat[:,1], label='Mathematica p=0 Rainbow')	
	#except:
	#	pass
				
else:
	if not os.path.exists("./data/mat_2nd_G_reducible"):
		os.system("math -script 2nd_G_Reducible.m")
	try:
		mat=np.loadtxt("data/mat_2nd_G_reducible")
		plt.plot(mat[:, 0], mat[:,1], label='Mathematica p=0 reducible')	
	except:
		pass
		
"""

nozero = np.nonzero(data2[:,1])
data2 = data2[nozero]
plt.errorbar(data2[:, 0], data2[:, 1], data2[:, 2], label='Sampled 2b')

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'second order a')
plt.legend()
plt.savefig('ana/plot_2a.pdf')
plt.clf()

#small plot
plt.errorbar(datasmall[:smend,0],datasmall[:smend,1], datasmall[:smend,2], label="Sampled")
try:
	plt.plot(mat[:smend, 0], mat[:smend, 1], label='Mathematica')
except:
	pass
plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G(\tau)$')
plt.legend()
plt.savefig('ana/plot_2a_small.pdf')
plt.clf()

#----------------- plot 2nd order b

#sample
data = np.loadtxt('data/second_orderb')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled 2b')

"""
if params['Self_Energy']:
	if not os.path.exists("./data/mat_2nd_g0se_crossed"):
		os.system("math -script 2nd_G0SE_Crossed.m")
	try:
		mat=np.loadtxt("data/mat_2nd_g0se_crossed")
		plt.plot(mat[:, 0], mat[:,1], label='Mathematica p=0 Crossed')	
	except:
		pass
else:
	if not os.path.exists("./data/mat_2nd_g_crossed"):
		os.system("math -script 2nd_G_Crossed.m")
	try:
		mat=np.loadtxt("data/mat_2nd_g_crossed")
		plt.plot(mat[:, 0], mat[:,1], label='Mathematica p=0 Crossed')	
	except:
		pass
"""		

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'second order b')
plt.legend()
plt.savefig('ana/plot_2b.pdf')
plt.clf()


	

#--------------- plot Ep vs ws

#data read in and sorting
Epvsws = np.loadtxt('data/Ep/Epvsws')
#Mathematica for 1st Order
if not os.path.exists("./data/Ep/mat_Ep1st"):
	os.system("math -script 1st_Ep.m")

Ep1st=np.loadtxt("data/Ep/mat_Ep1st")
for i in range(Epvsws[0,:].size - 1):
	Epvsws[:,i+1] += Ep1st[:,1]

#Ep vs Order
if params['SECumul']:
	Epvsord = np.loadtxt('data/Ep/Epvsord')
	for i in range(Epvsord[:,0].size):
		Epvsord[i, 1:] += Ep1st[:,1]

sort = np.argsort(-Epvsws[:,0])
Epvsws = Epvsws[sort]

#Pre Selection Estimator Data
Epmean = np.mean(Epvsws[:, 1:], axis=1)
Eperr = np.std(Epvsws[:, 1:], axis=1)/np.sqrt(len(Epvsws[0,:])-1)
Epdels = []
for Epit in range(Epmean.size):
	if (Eperr[Epit] > 0.25*abs(Epmean[Epit]) or (abs(Epmean[Epit]) - abs(params["Chemical_Potential"])) > abs(params["Chemical_Potential"])) and (abs(Epmean[Epit]) - abs(params["Chemical_Potential"])) > 1 :
		Epdels.append(Epit)
Epvsws=np.delete(Epvsws, Epdels, 0)
Epmean=np.delete(Epmean, Epdels)
Eperr=np.delete(Eperr, Epdels)

if params['SECumul']:	
	orders = Epvsord[:,0]
	Epvsordata = Epvsord[:,1:].T
	Epvsordata = Epvsordata[sort]
	Epvsordata = np.delete(Epvsordata, Epdels, 0)

#Plot Ep vs ws
# ws Gerade
def ws_gerade(ws):
	return ws
peter = np.vectorize(ws_gerade)

if Epmean.size == 0 :
	print('PreSelection failed!')

if Epmean.size > 1:
	plt.plot(-Epvsws[:, 0], peter(-Epvsws[:, 0]), 'r-', label=r'$\omega_{pol}$')
	plt.errorbar(-Epvsws[:, 0], Epmean, Eperr, label="Estimator")
	plt.xlabel(r'$-\omega_{pol}$')
	plt.ylabel(r'$-E_{pol}$')
	plt.legend()
	plt.savefig('ana/Epvsws.pdf')
	plt.clf()

	#finding root
	if params["Ep_Root"]:
		if params['SECumul']:
			Epa.Ep_root(Epvsws, Epvsordata, orders)
		else:
			Epa.Ep_root(Epvsws)
			

#--------------- plot Self Energy Estimator
#sample
data = np.loadtxt('data/SE/SE_1')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')

#Mathematica
if not os.path.exists("./data/SE/mat_1st_se"):
	os.system("math -script 1st_SE.m")
try:
	mat=np.loadtxt("data/SE/mat_1st_se")
	mat = mat[nozero]
	plt.plot(mat[:, 0], mat[:,1], label='Mathematica')	
except:
	pass

plt.yscale(SEscale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma^{1}(0, \tau)$')
plt.legend()
plt.savefig('ana/SE/SE_1.pdf')
plt.clf()

#small plot
plt.errorbar(data[:smend,0],data[:smend,1], data[:smend,2], label="Sampled")
try:
	plt.plot(mat[:smend, 0], mat[:smend, 1], label='Mathematica')
except:
	pass
plt.yscale(SEscale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma^{1}(0, \tau)$')
plt.legend()
plt.savefig('ana/SE/SE_1_small.pdf')
plt.clf()

#-------------------------------SE >1
#sample
data = np.loadtxt('data/SE/SE_>1')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')

#Mathematica
if params["Max_Order"]==2 or params["Total_Max_Order"]==2:
	if not os.path.exists("./data/SE/mat_2nd_se"):
		os.system("math -script 2nd_SE.m")
	try:
		mat=np.loadtxt("data/SE/mat_2nd_se")
		mat = mat[nozero]
		plt.plot(mat[:, 0], mat[:,1], label='Mathematica p=0')	
	except:
		pass

plt.yscale(SEscale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma^{>1}(0, \tau)$')
plt.legend()
plt.savefig('ana/SE/SE_>1.pdf')
plt.clf()

#small plot
plt.errorbar(data[:smend,0],data[:smend,1], data[:smend,2], label="Sampled")
plt.yscale(SEscale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma^{>1}(0, \tau)$')
plt.legend()
plt.savefig('ana/SE/SE_>1_small.pdf')
plt.clf()


#-------------------------------SE all
#sample
data = np.loadtxt('data/SE/SE_all')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')

plt.yscale(SEscale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma(0, \tau)$')
plt.legend()
plt.savefig('ana/SE/SE_all.pdf')
plt.clf()

#small plot
plt.errorbar(data[:smend,0],data[:smend,1], data[:smend,2], label="Sampled")
try:
	plt.plot(mat[:smend, 0], mat[:smend, 1]+data2[:smend,1], label='Mathematica')
except:
	pass
plt.yscale(SEscale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma(0, \tau)$')
plt.legend()
plt.savefig('ana/SE/SE_all_small.pdf')
plt.clf()

	
#-------------------------------------SE Transformation	
#Fitting G		
def Gtaugrwpse(tau, Z, Ek):
		return Z*np.exp(-(Ek-params["Chemical_Potential"])*tau)

if params["Self_Energy"]:
	if not os.path.exists("./data/mat_1st_g0se"):
		os.system("math -script 1st_G0SE.m")
	tause, setrans, taug0seiw, g0seiwtrans, taug0se, g0setrans = Trafos.green_trafo()
	plt.plot(taug0se, g0setrans, 'rx', label=r'$G_0 \Sigma$')
else:
	tause, setrans, taug0seiw, g0seiwtrans = Trafos.green_trafo()
	g0setrans = np.loadtxt('data/all_orders')
	plt.errorbar(g0setrans[:,0], g0setrans[:,1], g0setrans[:,2], label=r'G(\tau)')
	taug0se = g0setrans[:,0]
	g0setrans = g0setrans[:,1]
	
plt.plot(tause, setrans, 'b+', label=r'$\Sigma$')
plt.plot(taug0seiw, g0seiwtrans, 'g.', label=r'$G_0\Sigma$ Sampled')

#G0SE Fit
fitstart = -150
fitend = -10
Ep,Eperr, Z,_ = Epa.Ep_fit(fitstart,fitend, taug0se, g0setrans)
if Ep!=0:
	Fitdata = Epa.Gtaugrwp(taug0se,Z, Ep)
	nonirel = np.where(Fitdata> 1*10**(params['Irrelevance_Order']))
	plt.plot(taug0se[nonirel], Fitdata[nonirel], 'r',  label ='Ep = %g +- %.e' %(Ep, Eperr))
	
#SE Fit
fitstart = -60
fitend = -10
Ep,Eperr, Z,_ = Epa.Ep_fit(fitstart,fitend, tause, setrans)
if Ep!=0:
	Fitdata = Epa.Gtaugrwp(tause, Z,Ep)
	nonirel = np.where(Fitdata> 1*10**(params['Irrelevance_Order']))
	plt.plot(tause[nonirel], Fitdata[nonirel],'b',  label ='Ep = %g +- %.e' %(Ep, Eperr))
	
#G0SEiw Fit
fitstart = -200
fitend = -35
Ep,Eperr, Z,_ = Epa.Ep_fit(fitstart,fitend, taug0seiw, g0seiwtrans)
if Ep != 0:
	Fitdata = Epa.Gtaugrwp(taug0seiw,Z, Ep)
	nonirel = np.where(Fitdata> 10**(params['Irrelevance_Order']))
	plt.plot(taug0seiw[nonirel], Fitdata[nonirel], 'g', label ='Ep = %g +- %.e' %(Ep, Eperr))

plt.yscale('log')
plt.xlabel(r'$\tau$')
plt.ylabel(r'G($\tau$)')
plt.legend()
plt.savefig('ana/plot_Trans_all.pdf')
plt.clf()
	
	
	


