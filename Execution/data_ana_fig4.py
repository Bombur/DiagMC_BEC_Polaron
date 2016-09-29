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

# definitions
pi = math.pi
scale= 'log'
    
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
		  
silentremove('ana/plot_0.pdf')
silentremove('ana/plot_all.pdf')
silentremove('ana/plot_1.pdf')
silentremove('ana/plot_2a.pdf')
silentremove('ana/plot_2b.pdf')
silentremove('ana/Epvsws.pdf')
silentremove('ana/PolaronEnergy')
silentremove('ana/control/plot_q.pdf')
silentremove('ana/control/plot_all.pdf')
silentremove('ana/control/SEexpmut.pdf')
silentremove('ana/control/plot_first_Norm_End.pdf')
silentremove('ana/control/plot_last_Norm_End.pdf')



#--------------- plot Q Statistics
#sample
data = np.loadtxt('data/stats/qs_tot', skiprows=1)
norm = data[:,1].sum()
plt.plot(data[:, 0], data[:, 1]/norm*100)
plt.xlabel(r'Phonon Momentum $q$')
plt.ylabel(r'$P(q)$ [%]')
plt.xlim([0, params["Q_Cutoff"]])
plt.savefig('ana/control/plot_q.pdf')

plt.clf()

#--------------- plot Order Statistics
#sample
data = np.loadtxt('data/stats/os_tot', skiprows=1)
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.plot(data[:, 0], data[:, 1])
plt.xlabel(r'Order')
plt.ylabel(r'Counts')
plt.savefig('ana/control/plot_ord.pdf')

plt.clf()

#Loading Ws
data = np.loadtxt('data/Ep/Epvsws')
ws = data[:,0]

#SECUMUL
minmax= np.loadtxt('data/stats/minmax', skiprows = 3)
mimarows = minmax.size / 7
if params["SECumul"] and params["Self_Energy"] and mimarows >2:
	#Compare first Norm and End Diagram 
	data = np.loadtxt('data/secumul/Norm_first')
	nozero = np.nonzero(data[:,1])
	data = data[nozero]
	tmp = "%i" %(minmax[1,0])
	plt.plot(data[:, 0], data[:, 1], label=('Norm Diagram Order '+tmp))
	data2 = np.loadtxt('data/secumul/End_first')
	nozero = np.nonzero(data2[:,1])
	data2 = data2[nozero]	
	tmp = "%i" %(minmax[0,1])
	plt.plot(data2[:, 0], data2[:, 1], label='End Diagram Order ' +tmp)

	plt.yscale('log')
	plt.xlabel(r'$\tau$')
	plt.ylabel(r'Order '+tmp)
	#plt.xlim([0, taumax])
	#plt.ylim(auto=True)
	plt.legend()
	plt.savefig('ana/control/plot_first_Norm_End.pdf')
	plt.clf()
	
	#Compare last Norm and End Diagram 
	data = np.loadtxt('data/secumul/Norm_last')
	nozero = np.nonzero(data[:,1])
	data = data[nozero]	
	tmp = "%i" %(minmax[-1,0])
	plt.plot(data[:, 0], data[:, 1], label='Norm Diagram Order '+tmp)
	data2 = np.loadtxt('data/secumul/End_last')
	nozero = np.nonzero(data2[:,1])
	data2 = data2[nozero]	
	tmp = "%i" %(minmax[-2,1])
	plt.plot(data2[:, 0], data2[:, 1], label='End Diagram Order ' +tmp)

	plt.yscale('log')
	plt.xlabel(r'$\tau$')
	plt.ylabel(r'Order ' + tmp)
	#plt.xlim([0, taumax])
	#plt.ylim(auto=True)
	plt.legend()
	plt.savefig('ana/control/plot_last_Norm_End.pdf')
	plt.clf()
	
	#--------------- plot Normfac Statistics
	#sample
	preforder= minmax[:, 1]
	pref = np.empty_like(minmax[:,2])
	pref[0] = 1.
	for i in range(minmax[:,0].size-1):
		pref[i+1]= pref[i]*minmax[i,3]/minmax[i+1,2] 
	 
	plt.plot(preforder, pref)
	plt.yscale('log')
	plt.xlabel(r'Order')
	plt.ylabel(r'Norm Factor')
	plt.savefig('ana/control/plot_normfac.pdf')

	plt.clf()
	


#--------- plot Green-function

#sample
data = np.loadtxt('data/all_orders')
nozero = np.nonzero(data[:,1])
data = data[nozero]

fitFile = open("ana/PolaronEnergy", "w")
fitFile.write("#Polaron Energy \n")
fitFile.write("Ep \t Ep Error \t Z \t Z Error \n")

if params["Ep_Fit"]:
	#fitting data
	def Gtaugrwp(tau, Z, Ek):
		return Z*np.exp(-(Ek-params["Chemical_Potential"])*tau)
	
	nozero = np.nonzero(data[:,2])
	fitdata = data[nozero]
	
	try:
		popt, pcov = optimize.curve_fit(Gtaugrwp, fitdata[params["Fit_Bin"]:-5,0], fitdata[params["Fit_Bin"]:-5,1], p0= [1,params["Chemical_Potential"]+1], sigma= fitdata[params["Fit_Bin"]:-5,2] )
		perr = np.sqrt(np.diag(pcov)) 
		if perr[1] < 0.25*abs(popt[1]) or (abs(popt[1] - params["Chemical_Potential"])) < abs(params["Chemical_Potential"]):
			print("Accepted Fit!")
			fitFile.write("%.3e \t %.3e \t %.3e \t %.3e\n" %(popt[1], perr[1], popt[0], perr[0]))
			fitcurve = np.vectorize(Gtaugrwp)
			plt.plot(data[:, 0], Gtaugrwp(data[:,0], popt[0], popt[1]), label="Fit alpha="+str(params["Alpha"])+" E=%.2f" %(popt[1]))
		
	except:
		print("Fit did not work!")
fitFile.close()

plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')			
plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma (\tau)$ or $G(\tau)$')
plt.legend()
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
plt.xlabel(r'$\tau$')
plt.ylabel(r'Fake-function')
plt.legend()
plt.savefig('ana/plot_0.pdf')
plt.clf()



#--------------- plot 1st order
#sample
data = np.loadtxt('data/first_order')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')

# Compare
if params["Mat_G0SEG0"]:
	if not os.path.exists("./data/mat_1st_g0seg0"):
		os.system("math -script BEC_1st_G0SEG0.m")
	try:
		peter=np.loadtxt("data/mat_1st_g0seg0")
		plt.plot(peter[:, 0], peter[:,1], label='Mathematica')	
	except:
		os.system("math -script BEC_1st_G0SEG0.m")
		peter=np.loadtxt("data/mat_1st_g0seg0")
		plt.plot(peter[:, 0], peter[:,1], label='Mathematica')	
	

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_0 \Sigma G_0 (0, \tau)$')
plt.legend()
plt.savefig('ana/plot_1.pdf')
plt.clf()

#--------------- plot 2nd order a
#sample
data = np.loadtxt('data/second_ordera')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled 2a')
data2 = np.loadtxt('data/second_orderb')
nozero = np.nonzero(data2[:,1])
data2 = data2[nozero]
plt.errorbar(data2[:, 0], data2[:, 1], data2[:, 2], label='Sampled 2b')

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'second order a')
plt.xlim([0, taus[2]])
plt.ylim(auto=True)
plt.legend()
plt.savefig('ana/plot_2a.pdf')
plt.clf()

#----------------- plot 2nd order b

#sample
data = np.loadtxt('data/second_orderb')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled 2b')

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'second order b')
plt.xlim([0, taus[2]])
plt.ylim(auto=True)
plt.legend()
plt.savefig('ana/plot_2b.pdf')
plt.clf()


	

#--------------- plot Ep vs ws

#data read in and sorting
#First Order in Mathematica
if not os.path.exists("./data/Ep/mat_1st_Epvsws"):
	os.system("math -script BEC_1st_SEomega.m")
Ep1st = np.loadtxt('data/Ep/mat_1st_Epvsws')
sort = np.argsort(-Ep1st[:,0])
Ep1st = Ep1st[sort]
#Data
Epvsws = np.loadtxt('data/Ep/Epvsws')
Epvsws = Epvsws[sort]
#Ep vs Order
Epvsord = np.loadtxt('data/Ep/Epvsord')
orders = Epvsord[:,0]
Epvsordata = Epvsord[:,1:].T
Epvsordata = Epvsordata[sort]

#Pre Selection Estimator Data
Epmean = np.mean(Epvsws[:, 1:], axis=1)
Eperr = np.std(Epvsws[:, 1:], axis=1)/np.sqrt(len(Epvsws[0,:])-1)
Epdels = []
for Epit in range(Epmean.size):
	if (Eperr[Epit] > 0.25*abs(Epmean[Epit]) or (abs(Epmean[Epit]) - abs(params["Chemical_Potential"])) > abs(params["Chemical_Potential"])) and (abs(Epmean[Epit]) - abs(params["Chemical_Potential"])) > 1 :
		Epdels.append(Epit)
Epvsws=np.delete(Epvsws, Epdels, 0)
Ep1st=np.delete(Ep1st, Epdels, 0)
Epvsordata = np.delete(Epvsordata, Epdels, 0)
Epmean=np.delete(Epmean, Epdels)
Eperr=np.delete(Eperr, Epdels)


#Plot Ep vs ws
# ws Gerade
def ws_gerade(ws):
	return ws
peter = np.vectorize(ws_gerade)

if Epmean.size == 0 :
	print('PreSelection failed!')

if Epmean.size > 1:
	plt.plot(-Ep1st[:, 0], peter(-Ep1st[:, 0]), 'r-', label=r'$\omega_{pol}$')
	plt.plot(-Ep1st[:, 0], -Ep1st[:, 1], 'go', label="1st Order")
	plt.errorbar(-Epvsws[:, 0], Epmean, Eperr, label="Estimator")
	plt.xlabel(r'$-\omega_{pol}$')
	plt.ylabel(r'$-E_{pol}$')
	plt.legend()
	plt.savefig('ana/Epvsws.pdf')
	plt.clf()

	#finding root
	if params["Ep_Root"]:
		Eproots = []
		x = - Epvsws[:,0] 
		for col in Epvsws[:,1:].T: 
			y = col + Epvsws[:,0]
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
			
			#Check Epvs Ord to control
			if abs(col[aidx]-Epvsordata[aidx, -1]) > col[aidx]:
				continue
			if abs(Epvsordata[aidx, -1] - Epvsordata[aidx, -2]) > 0.001*Epvsordata[aidx, -1]:
				print('Estimator did not converge!')
				continue
			plt.plot(orders, Epvsordata[aidx], label  = r'$\omega_{lowl} = -%g$' %(a))
			plt.plot(orders, Epvsordata[bidx], label  = r'$\omega_{uppl} = -%g$' %(b))
			
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
		
		plt.xlabel('Order')
		plt.ylabel(r'$-E_{pol}$')
		plt.legend()
		plt.savefig('ana/control/Epvsord.pdf')
		plt.clf()			
		
		if len(Eproots) == 1:
			with open('ana/PolaronEnergy', 'r+') as fitin:
				fitin.seek(0,2)
				fitin.write("%.3e \n" %(Eproots[0]))
		elif len(Eproots) > 1:
			with open('ana/PolaronEnergy', 'r+') as fitin:
				fitin.seek(0,2)
				fitin.write("%.3e \t %.3e \n" %(np.mean(Eproots), np.std(Eproots)/np.sqrt(len(Eproots)-1)))


#--------------- plot Ep Integrand Statistics

#data = np.loadtxt('data/Ep/Eptest')
#iti=0
#for col in data.T[1:]:
#	plt.plot(data[:, 0], col, label=r'$\omega_{Pol} = %g$' %(ws[iti]))
#	iti += 1
#peter = np.loadtxt(peter_loc+ '/mat_FP_1st_seexpmut_w2')
#plt.plot(peter[:,0], peter[:,1], label="Mat") 
#plt.yscale('log')
#plt.xlabel(r'$\tau$')
#plt.ylabel(r'$\Sigma(\mathbf{p}=0, \tau) e^{(\omega-\mu) \tau} $')
#plt.xlim([0, 1.5])
#plt.ylim([1e-2,1000])
#plt.legend()
#plt.savefig('ana/control/SEexpmut.pdf')

#plt.clf()
	
	
	
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
		plt.xlim([0, taumax])
		plt.legend()
		plt.savefig('ana/control/plot_all.pdf')
		plt.clf()

	



