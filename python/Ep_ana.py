#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import json
import os, errno, sys
from scipy import optimize
from scipy.interpolate import interp1d
from multiprocessing import Pool

#define constants
pi = np.pi

# import parameters
with open('DiagMC_BEC.json') as params_file:
    params = json.load(params_file)
p = params['Momentum']
Eprenorm = params['Alpha']/np.sqrt(2)/pi *(1+ 1/params['Impurity_Mass'])*params['Q_Cutoff']
print(Eprenorm)

#------------------------Fit Function
def Gtaugrwp(tau, Z, Ek):
	return Z*np.exp(-(Ek+ p**2/2 -params["Chemical_Potential"])*tau)

def Ep_fit(fitstart, fitend, tau, data, err = np.array([])):
	"""
	Fit data in tau to exponential function
	
    Args:
    	fitstart and fitend: indexes to fit inbetween 
        tau: Numpy array 
        data: Numpy array
        err: Numpy array with error (optional)
        
    Returns:
        4 double Ep, Error of Ep, Z and Error of Z
        also to PolaronEnergy file
	"""
	try:
		if err.size != 0:
			popt, pcov = optimize.curve_fit(Gtaugrwp, tau[fitstart:fitend], data[fitstart:fitend], p0= [1, params["Chemical_Potential"]], sigma=err[fitstart:fitend] )
		else:
			popt, pcov = optimize.curve_fit(Gtaugrwp, tau[fitstart:fitend], data[fitstart:fitend], p0= [1, params["Chemical_Potential"]])
		
		perr = np.sqrt(np.diag(pcov)) 
		
		if perr[1] < 0.25*abs(popt[1]) or (abs(popt[1] - params["Chemical_Potential"])) < abs(params["Chemical_Potential"]):
			print("Accepted Fit!")
			if params['BEC_Polaron']:
				with open('ana/PolaronEnergy', 'r+') as fitin:
					fitin.seek(0,2)
					fitin.write("%.3e \t %.3e \t %.3e \t %.3e\n" %(Eprenorm+popt[1], perr[1], popt[0], perr[0]))
				return popt[1], perr[1], popt[0], perr[0]
			else:
				with open('ana/PolaronEnergy', 'r+') as fitin:
					fitin.seek(0,2)
					fitin.write("%.3e \t %.3e \t %.3e \t %.3e\n" %(popt[1], perr[1], popt[0], perr[0]))
				return popt[1], perr[1], popt[0], perr[0]
		else:
			print("Fit not accepted!")
			return 0,0,0,0
	except:
		print("Fit did not work!")
		return 0,0,0,0



#---------------------Find_Root
def Ep_root(Epvsws, Epvsord=0, orders=0):
	"""
	Find intersection of the Estimator data and the preselected omegas

    Args:
        Epvsws: Numpy array 1.col: preselected omegas
        					rightcols: Estimator data
        Epvsord: Numpy array in case of cumulative sampling
        orders: Array of order Steps
        
    Returns:
        Ep written to PolaronEnergy file
    """
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
		if params['SECumul'] and Epvsord:
			plt.plot(orders, Epvsord[aidx], label  = r'$\omega_{lowl} = -%g$' %(a))
			plt.plot(orders, Epvsord[bidx], label  = r'$\omega_{uppl} = -%g$' %(b))
			if abs(col[aidx]-Epvsord[aidx, -1]) > col[aidx]:
				continue
			if abs(Epvsord[aidx, -1] - Epvsord[aidx, -2]) > 0.00001*Epvsord[aidx, -1]:
				print('Estimator did not converge!')
				continue
		
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
	plt.legend(loc=4)
	plt.savefig('ana/Ep/Epvsord.pdf')
	plt.clf()
		
	if len(Eproots) == 1:
		with open('ana/PolaronEnergy', 'r+') as fitin:
			fitin.seek(0,2)
			if params['BEC_Polaron']:
				fitin.write("%.5e \n" %(-Eproots[0]+Eprenorm))
			else:
				fitin.write("%.5e \n" %(-Eproots[0]))
	elif len(Eproots) > 1:
		with open('ana/PolaronEnergy', 'r+') as fitin:
			fitin.seek(0,2)
			if params['BEC_Polaron']:
				fitin.write("%.6e \t %.4e \n" %(-np.mean(Eproots)+Eprenorm, np.std(Eproots)/np.sqrt(len(Eproots)-1)))
			else:
				fitin.write("%.6e \t %.4e \n" %(-np.mean(Eproots), np.std(Eproots)/np.sqrt(len(Eproots)-1)))
				


#Statistics Analysis
def Stat_ana():
	"""
	#---------------Test histo
	data= np.loadtxt('data/stats_sswap/testhisto')
	plt.plot(data[:,0]/data[:,0].sum(), label='left')
	plt.plot(data[:,1]/data[:,1].sum(), label='right')
	data= np.loadtxt('data/stats/testhisto')
	plt.plot(data[:,2]/data[:,2].sum(), label='meas')
	plt.xlabel(r'$\tau$')
	plt.xlim([0,50])
	plt.ylabel(r'$P(\tau)$ [%]')
	plt.legend()
	plt.yscale('log')
	plt.savefig('ana/control/plot_testhisto.pdf')
	plt.clf()
	"""


	#--------------- plot Q Statistics
	#all
	data = np.loadtxt('data/stats/qs_tot', skiprows=1)
	data1=data
	nozero = np.nonzero(data[:,1])
	data= data[nozero] 
	norm = data[:,1].sum()
	data[:,1] /=norm
	cumul= np.zeros_like(data[:,1])
	cumul[0]= data[0,1]
	for qprop in range(data[:,1].size -1): 
		cumul[qprop +1] = cumul[qprop] + data[qprop+1,1]
	plt.plot(data[:, 0], cumul)
	plt.xlabel(r'Phonon Momentum $q$')
	plt.ylabel(r'$P(q)$ [%]')
	plt.savefig('ana/control/plot_q.pdf')
	plt.clf()
	
	#1
	nozero = np.nonzero(data1[:,2])
	data1= data1[nozero] 
	data1[:,2] /= data1[:,2].sum()*100
	plt.plot(data1[:, 0], data1[:,2])
	plt.xlabel(r'Phonon Momentum $q$')
	plt.ylabel(r'$P(q)$ [%]')
	plt.savefig('ana/control/plot_q_1.pdf')
	plt.clf()
	
	
	
	#---------------- plot tau Statistics
	data = np.loadtxt('data/stats/taus_tot', skiprows=1)
	nozero = np.nonzero(data[:,1])
	dataall= data[nozero] 
	norm = dataall[:,1].sum()
	dataall[:,1] = dataall[:, 1]/norm *100
	plt.plot(dataall[:, 0], dataall[:,1], label="All")
	
	#1
	nozero = np.nonzero(data[:,2])
	data1 = data[nozero]
	data1[:,2] = data1[:, 2]/data1[:,2].sum() *100
	plt.plot(data1[:, 0], data1[:,2], label="1")
	
	#crossed
	nozero = np.nonzero(data[:,3])
	datac = data[nozero]
	datac[:,3] = datac[:, 3]/datac[:,3].sum() *100
	plt.plot(datac[:, 0], datac[:,3],label="Crossed")
	
	#Rainbow
	nozero = np.nonzero(data[:,4])
	datar = data[nozero]
	datar[:,4] = datar[:, 4]/datar[:,4].sum() *100
	plt.plot(datar[:, 0], datar[:,4], label="Rainbow")
	plt.yscale('log')
	plt.legend()
	plt.xlabel(r'$\tau$')
	plt.ylabel(r'$P(\tau)$ [%]')
	plt.savefig('ana/control/plot_tau.pdf')
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
	
	
	#--------------------------SECUMUL
	if params["SECumul"]:
		minmax= np.loadtxt('data/stats/minmax', skiprows = 3, usecols=(0,1,2,3))
		mimarows = minmax.size/4
	
		if mimarows>2:
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
			plt.legend()
			plt.savefig('ana/control/plot_last_Norm_End.pdf')
			plt.clf()
			
			#--------------- plot Normfac Statistics
			#sample
			preforder = minmax[:, 0]
			pref = np.empty(mimarows)
			pref[0] = 1.
			for i in range(minmax[:,0].size-1):
				pref[i+1]= pref[i]*minmax[i,3]/minmax[i+1,2] 
	
			plt.plot(preforder, pref)
			plt.yscale('log')
			plt.xlabel(r'Order')
			plt.ylabel(r'Norm Factor')
			plt.savefig('ana/control/plot_normfac.pdf')
			plt.clf()
		
	
#----------------------run Mathematica Scripts
def mat_scripts():
	pool= Pool(processes = params['NCores'])
	matalways = ['1st_Ep.m','1st_SE.m','1st_G0SEG0.m',]
	matgsamp = []
	matsesamp = ['1st_G0SE.m']
	if params['Second_Order_Mathematica']:
		matgsamp.append('2nd_G_Reducible.m')
		matgsamp.append('2nd_G_Crossed.m')
		matsesamp.append('2nd_G0SE_Crossed.m')
		matsesamp.append('2nd_G0SE_Rainbow.m')
		matsesamp.append('2nd_G0SE_Complete.m')
	if params['Mat_Scripts']:
		for i in matalways:
			command='math -script ' + i
			pool.apply_async(os.system, (command,))
		if params['Self_Energy']:
			for i in matsesamp:
				command='math -script ' + i
				pool.apply_async(os.system, (command,))
		else:
			pass
			for i in matgsamp:
				command='math -script ' + i
				pool.apply_async(os.system, (command,))
				
		pool.close()
		pool.join()

