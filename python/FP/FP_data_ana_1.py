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

# import parameters
with open('DiagMC_BEC.json') as params_file:
    params = json.load(params_file)

taus = params['Taus']
taumax = taus[-1]
bins = params['Bins']
peter_loc = params["Peter_Path"]
p  = 0.3
#Loading Ws
data = np.loadtxt('data_FP_ord1/Ep/Epvsws')
ws = data[:,0]

# definitions
pi = math.pi
scale= 'linear'
    
#helper function
def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

# cleanup
def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured
		  


#---------------- plot Green-function

#sample
data = np.loadtxt('data_FP_ord1/all_orders')
nozero = np.nonzero(data[:,1])
data = data[nozero]

fitFile = open("ana_FP_ord1/PolaronEnergy", "w")
fitFile.write("#Polaron Energy \n")
fitFile.write("Ep \t Ep Error \t Z \t Z Error \n")

if params["Ep_Fit"]:
	#fitting data
	def Gtaugrwp(tau, Z, Ek):
		return Z*np.exp(-(Ek +(p*p)/2-params["Chemical_Potential"])*tau)
	
	nozero = np.nonzero(data[:,2])
	fitdata = data[nozero]
	
	#params['Fit_Bin'] = 50
	fitend = -10
	try:
		popt, pcov = optimize.curve_fit(Gtaugrwp, fitdata[params["Fit_Bin"]:fitend,0], fitdata[params["Fit_Bin"]:fitend,1], p0= [2,params["Chemical_Potential"]], sigma= fitdata[params["Fit_Bin"]:fitend,2] )
		perr = np.sqrt(np.diag(pcov)) 
		if perr[1] < 0.25*abs(popt[1]) or (abs(popt[1] - params["Chemical_Potential"])) < abs(params["Chemical_Potential"]):
			print("Accepted Fit!")
			fitFile.write("%.3e \t %.3e \t %.3e \t %.3e\n" %(popt[1], perr[1], popt[0], perr[0]))
			plt.plot(data[:, 0], Gtaugrwp(data[:,0], popt[0], popt[1]),'b', label="Ep=%.2f +- %.e " %(popt[1], perr[1]))
		else:
			print("Fit not accepted!")
	except:
		print("Fit did not work!")
fitFile.close()

plt.errorbar(data[:, 0], data[:, 1], data[:, 2], fmt = 'bx', label='Sampled')		
peter = np.loadtxt(peter_loc + '/FP_MO1')

try:
	popt, pcov = optimize.curve_fit(Gtaugrwp, peter[params["Fit_Bin"]:,0], peter[params["Fit_Bin"]:,1], p0= [2,params["Chemical_Potential"]])
	perr = np.sqrt(np.diag(pcov)) 
	
except:
	print("Fit did not work!")
plt.plot(peter[:,0], peter[:,1], 'rx', label ='Peter p=0.3')	
plt.plot(peter[:,0], Gtaugrwp(peter[:,0], popt[0], popt[1]), 'r', label ='Peter Ep = %.2f +- %.e' %(popt[1], perr[1]))


plt.yscale('log')
plt.xlim([0,10])
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G(\tau$')
plt.legend(loc = 3)
plt.savefig('ana_FP_ord1/plot_all.pdf')

plt.clf()


#--------------- plot Ep vs ws

#data read in and sorting
Epvsws = np.loadtxt('data_FP_ord1/Ep/Epvsws')
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


#Plot Ep vs ws
# ws Gerade
def ws_gerade(ws):
	return ws + p*p/2
peter = np.vectorize(ws_gerade)

if Epmean.size == 0 :
	print('PreSelection failed!')

if Epmean.size > 1:
	plt.plot(-Epvsws[:, 0], peter(-Epvsws[:, 0]), 'r-', label=r'$\omega_{pol}$')
	plt.errorbar(-Epvsws[:, 0], Epmean, Eperr, label="Estimator")
	plt.xlabel(r'$-\omega_{pol}$')
	plt.ylabel(r'$-E_{pol}$')
	plt.legend()
	plt.savefig('ana_FP_ord1/Epvsws.pdf')
	plt.clf()

	#finding root
	if params["Ep_Root"]:
		Eproots = []
		x = - Epvsws[:,0] 
		for col in Epvsws[:,1:].T: 
			y = col + Epvsws[:,0] - p*p/2
			func = interp1d(x, y)
			if not (np.amin(y) < 0 and np.amax(y) > 0) :
				print('Ep out of range to find root!')
				continue
			ysmaller = np.where(y < 0)
			ysmaller = ysmaller[0]
			ybigger = np.where(y > 0)
			ybigger = ybigger[0]
			aidx = ysmaller[np.argmax(y[ysmaller])]
			bidx = ybigger[np.argmin(y[ybigger])]
			if (aidx - bidx) != 1:
				continue
			a = x[ysmaller[np.argmax(y[ysmaller])]]
			b = x[ybigger[np.argmin(y[ybigger])]]
			
			try:
				Eproots.append(optimize.brentq(func, a, b))
				print("brentq() converged!")
			except (ValueError, RuntimeError) as e:
				try:
					Eproots.append(optimize.brenth(func, a, b))
					print("brenth() converged!")
				except (ValueError, RuntimeError) as e:
					try:
						Eproots.append(optimize.ridder(func, a, b))
						print("ridder() converged!")
					except (ValueError, RuntimeError) as e:
						try:
							Eproots.append(optimize.bisect(func, a, b))
							print("bisect() converged!")
						except (ValueError, RuntimeError) as e:
							try:
								Eproots.append(optimize.newton(func, x[1]))
								print("newton() converged!")
							except (ValueError, RuntimeError) as e:
								print("No Estimator Result!")
								continue
		

		if len(Eproots) == 1:
			with open('ana_FP_ord1/PolaronEnergy', 'r+') as fitin:
				fitin.seek(0,2)
				fitin.write("%.3e \n" %(Eproots[0]))
		elif len(Eproots) > 1:
			with open('ana_FP_ord1/PolaronEnergy', 'r+') as fitin:
				fitin.seek(0,2)
				fitin.write("%.3e \t %.3e \n" %(np.mean(Eproots), np.std(Eproots)/np.sqrt(len(Eproots)-1)))



