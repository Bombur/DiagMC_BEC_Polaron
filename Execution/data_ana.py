#!/usr/bin/env python3
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, errno, sys
import math
import json
from scipy.optimize import curve_fit


# import parameters
with open('DiagMC_BEC.json') as params_file:
    params = json.load(params_file)

# data for the comparison    
peter_loc = params["Peter_Path"]
if params["Froehlich_Polaron"]:
	if params["Self_Energy" ] and params["Max_Order"] == 6:
		compare = ["FP_1","FP_2","FP_SE_6ord"]
	else:
		compare = ["FP_1","FP_2","FP_all"]
	scale = 'linear'
else:
	compare = ["mat_1st_g0seg0", "mat_2nd_g0se","fig3total"]
	scale= 'log'
    
# cleanup
def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured
		  
silentremove('ana/plot_0.pdf')
silentremove('ana/plot_all.pdf')
silentremove('ana/plot_1.pdf')
silentremove('ana/plot_2a.pdf')
silentremove('ana/plot_2b.pdf')
silentremove('ana/control/plot_first_Norm_End.pdf')
silentremove('ana/control/plot_last_Norm_End.pdf')


# definitions
pi = math.pi

# plot fake-function
#sample
fake = np.loadtxt('data/zero_order')
plt.errorbar(fake[:, 0], fake[:, 1], fake[:, 2], label='Sampled')

# analytical
def fake_integrand(tau):
  disp = np.power(params['Momentum'], 2)/2. - params['Chemical_Potential']
  return np.exp(- disp * tau)

fake_integrand = np.vectorize(fake_integrand)
plt.plot(fake[:, 0], fake_integrand(fake[:, 0]), label='Analytical result')
plt.xlabel(r'$\tau$')
plt.ylabel(r'Fake-function')
plt.xlim([0, params['Tau_max']])
plt.legend()
plt.savefig('ana/plot_0.pdf')

plt.clf()

# plot Green-function
#def mucor(tau):
#  return np.exp(params['Chemical_Potential'] * tau)
#mucor = np.vectorize(mucor)

#sample
data = np.loadtxt('data/all_orders')
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')

# Peter
peter = np.loadtxt(peter_loc+'/'+compare[2])
plt.plot(peter[:, 0], peter[:, 1], label=compare[2])

#fitting data
def Gtaugrwp(tau, Z, Ek):
	return Z*np.exp(-(Ek-params["Chemical_Potential"])*tau)

#write fit file
#popt, pcov = curve_fit(Gtaugrwp, data[50:,0], data[50:,1])
#perr = np.sqrt(np.diag(pcov)) 
#fitFile = open("ana/fitparams", "w")
#fitFile.write("#Fit Parameters \n")
#fitFile.write("Polaron Energy \n")
#fitFile.write("E = "+str(popt[1])+" +- "+ str(perr[1])+"\n \n")
#fitFile.write("Z = "+str(popt[0])+" +- "+ str(perr[0]))
#fitFile.close()

#fitcurve = np.vectorize(Gtaugrwp)
#plt.plot(data[:, 0], Gtaugrwp(data[:,0], popt[0], popt[1]), label="Fit alpha="+str(params["Coupling_Strength"])+" E=%.2f" %(popt[1]))
plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma (\tau) e^{\mu \tau}$')
plt.xlim([0, params['Tau_max']])
plt.legend()
plt.savefig('ana/plot_all.pdf')

plt.clf()

# plot 1st order
#sample
data = np.loadtxt('data/first_order')
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')

# Compare
peter = np.loadtxt(peter_loc+ '/'+compare[0])
plt.plot(peter[:, 0], peter[:, 1], label=compare[0])
plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_0 \Sigma G_0 (0, \tau)$')
plt.xlim([0, params['Tau_max']])
plt.legend()
plt.savefig('ana/plot_1.pdf')

plt.clf()

# plot 2nd order a
#sample
data = np.loadtxt('data/second_ordera')
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled 2a')
data2 = np.loadtxt('data/second_orderb')
plt.errorbar(data2[:, 0], data2[:, 1], data2[:, 2], label='Sampled 2b')

# Compare
peter = np.loadtxt(peter_loc+ '/'+ compare[1])
plt.plot(peter[:, 0], peter[:, 1], label=compare[1])
plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'second order a')
plt.xlim([0, params['Tau_max']])
plt.legend()
plt.savefig('ana/plot_2a.pdf')

plt.clf()

# plot 2nd order b
#sample
data = np.loadtxt('data/second_orderb')
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled 2b')

# Compare
plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'second order b')
plt.xlim([0, params['Tau_max']])
plt.legend()
plt.savefig('ana/plot_2b.pdf')
plt.clf()


#SECUMUL
if params["SECumul"] and params["Self_Energy"]:
	minmax= np.loadtxt('data/stats/minmax', skiprows = 3)
	
	#Compare first Norm and End Diagram 
	data = np.loadtxt('data/secumul/Norm_first')
	tmp = "%i" %(minmax[1,0])
	plt.plot(data[:, 0], data[:, 1], label=('Norm Diagram Order '+tmp))
	data2 = np.loadtxt('data/secumul/End_first')
	tmp = "%i" %(minmax[0,1])
	plt.plot(data2[:, 0], data2[:, 1], label='End Diagram Order ' +tmp)

	plt.yscale(scale)
	plt.xlabel(r'$\tau$')
	plt.ylabel(r'Order '+tmp)
	plt.xlim([0, params['Tau_max']])
	plt.legend()
	plt.savefig('ana/control/plot_first_Norm_End.pdf')
	plt.clf()
	
	#Compare last Norm and End Diagram 
	data = np.loadtxt('data/secumul/Norm_last')
	tmp = "%i" %(minmax[-1,0])
	plt.plot(data[:, 0], data[:, 1], label='Norm Diagram Order '+tmp)
	data2 = np.loadtxt('data/secumul/End_last')
	tmp = "%i" %(minmax[-2,1])
	plt.plot(data2[:, 0], data2[:, 1], label='End Diagram Order ' +tmp)

	plt.yscale(scale)
	plt.xlabel(r'$\tau$')
	plt.ylabel(r'Order ' + tmp)
	plt.xlim([0, params['Tau_max']])
	plt.legend()
	plt.savefig('ana/control/plot_last_Norm_End.pdf')
	plt.clf()
	
	#FP Check of Selfeneergy
if params["Self_Energy"]:
	if params["Froehlich_Polaron"]:
		#silentremove('ana/plot_all.pdf')
		
		if params["FOG0SE"]:
			os.system("math -script FP_g0SE_Transform.m")
		else:
			os.system("math -script FP_g0SEg0_Transform.m")
		
		#new all orders Green Function
		data = np.loadtxt('data/FP_control/mat_FP_all_transform')
		plt.plot(data[:, 0], data[:, 1], label='FP of SE Sampling Simpson')
		peter = np.loadtxt(peter_loc+'/'+compare[2])
		plt.plot(peter[:, 0], peter[:, 1], label=compare[2])
	
		plt.yscale(scale)
		plt.xlabel(r'$\tau$')
		plt.ylabel(r'G($\tau$)')
		plt.xlim([0, params['Tau_max']])
		plt.legend()
		plt.savefig('ana/control/plot_all.pdf')
		plt.clf()

	



