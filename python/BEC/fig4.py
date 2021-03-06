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
peter_loc = params["Peter_Path"] + '/'

# definitions
pi = math.pi
scale= 'log'

#helper function
def find_nearest(array,value):
    return (np.abs(array-value)).argmin()
    

#-------------------------------------Fitting
def Gtaugrwp(mu, tau, Z, Ek):
	return Z*np.exp(-(Ek-mu)*tau)
	

count=0
for dirpath, dirnames, filenames in os.walk(os.getcwd()):
	if "all_orders" in filenames:
		count +=1
fig4fit = np.empty([count, 3])
fig4dels = []

#control Plots for the fitting
alphalist = params["Alpha_List"]
alphalist = np.asarray(alphalist)
fig =[]
ax =[]
for j in range(len(alphalist)) :
	figtmp, axtmp = plt.subplots()
	fig.append(figtmp)
	ax.append(axtmp)

it=0
for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/data/all_orders"):
		os.chdir(folder)
		#Import Parameters of each run
		with open('DiagMC_BEC.json') as parit_file:
			parit = json.load(parit_file)
    	  						
    	#Fitting	  				
		def fitfunc(tau, Z, Ek):
			return Gtaugrwp(parit["Chemical_Potential"], tau, Z, Ek)
			
		data = np.loadtxt('data/all_orders')
		nozero = np.nonzero(data[:,1])
		data = data[nozero]
		nozero = np.nonzero(data[:,2])
		data = data[nozero]
		
		#parit["Fit_Bin"]+=10		

		try:
			popt, pcov = optimize.curve_fit(fitfunc, data[parit["Fit_Bin"]:,0], data[parit["Fit_Bin"]:,1], p0= [1,parit["Chemical_Potential"]+1], sigma= data[parit["Fit_Bin"]:,2] )
			perr = np.sqrt(np.diag(pcov)) 
			if perr[1] > 0.25*abs(popt[1]) or (abs(popt[1] - parit["Chemical_Potential"])) > abs(parit["Chemical_Potential"]):
				fig4dels.append(it)
				os.chdir(os.pardir)
				it += 1
				continue
			print("Accepted Fit!")			 		
		except:
			fig4dels.append(it)
			os.chdir(os.pardir)
			it += 1
			continue
			
		#Fit Data Transfer
		fig4fit[it, 0] = parit["Alpha"] 
		fig4fit[it, 1] = popt[1]
		fig4fit[it, 2] = perr[1]
		
		#Fit Control Plots
		idx = find_nearest(alphalist, parit["Alpha"])
		ax[idx].errorbar(data[:, 0], data[:, 1], data[:, 2], label='mu =  %.2e' %(parit["Chemical_Potential"]) )
		ax[idx].plot(data[:,0], fitfunc(data[:,0], popt[0], popt[1]) , label="E = %.2e" %(popt[1]))
		ax[idx].set_title("Alpha = " +str(parit["Alpha"]))
		os.chdir(os.pardir)
		it+=1


fig4fit=np.delete(fig4fit, fig4dels, 0)

for pit in range(len(ax)):
	ax[pit].set_yscale('log')
	ax[pit].set_xlabel(r'$\tau$')
	ax[pit].set_ylabel(r'$\Sigma (\mathbf{p} =0, \tau)$')
	ax[pit].legend(loc=2)
	fig[pit].savefig('all_%.1f.pdf' %(alphalist[pit]))
	ax[pit].cla()
	


	
#-----------------------------------------Estimator	

fig4esti = np.empty([len(alphalist), 3])
Eprootslist = []
Epconvord = []

it=0
for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/data/Ep/Epvsws"):
		os.chdir(folder)

		if not os.path.exists("./data/Ep/mat_1st_Epvsws"):
			os.system('math -script BEC_1st_SEomega.m')

		#Import Parameters of each run
		with open('DiagMC_BEC.json') as parit_file:
			parit = json.load(parit_file)
		
		#Pre Selection Estimator Data
		Epvsws = np.loadtxt('data/Ep/Epvsws')
		Ep1st = np.loadtxt('data/Ep/mat_1st_Epvsws')
		Epvsord = np.loadtxt('data/Ep/Epvsord')
		if Epvsord.ndim != 2 :
			os.chdir(os.pardir)
			it+=1
			continue
		orders = Epvsord[:,0]
		Epvsordata = Epvsord[:,1:].T
		sort = np.argsort(-Ep1st[:,0])
		Ep1st = Ep1st[sort]
		Epvsws = Epvsws[sort]
		Epvsordata = Epvsordata[sort]

		Epmean = np.mean(Epvsws[:, 1:], axis=1)
		Eperr = np.std(Epvsws[:, 1:], axis=1)/np.sqrt(len(Epvsws[0,:])-2)
		Epdels = []
		for Epit in range(Epmean.size):
			if Eperr[Epit] > 0.25*abs(Epmean[Epit]) or ((abs(Epmean[Epit]) - abs(parit["Chemical_Potential"])) > abs(parit["Chemical_Potential"]) and (abs(Epmean[Epit]) - abs(params["Chemical_Potential"])) > 1 ):
				Epdels.append(Epit)
		Epvsws=np.delete(Epvsws, Epdels, 0)
		Ep1st=np.delete(Ep1st, Epdels, 0)
		Epvsordata = np.delete(Epvsordata, Epdels, 0)

		rows, cols = Epvsws.shape
		if rows < 2 :
			os.chdir(os.pardir)
			it+=1
			continue
		
		#finding root
		Eproots = []
		x = -Epvsws[:,0]
		for col in Epvsws[:,1:].T:
			y = col+Epvsws[:,0]
			func = interp1d(x, y)
			if not (np.amin(y) < 0 and np.amax(y) > 0) :
				continue
			ysmaller = np.where(y < 0)
			ysmaller = ysmaller[0]
			ybigger = np.where(y > 0)
			ybigger = ybigger[0]
			aidx = ysmaller[np.argmax(y[ysmaller])]
			bidx = ybigger[np.argmin(y[ybigger])]
			if (aidx - bidx) != 1:
				continue
			a = x[aidx]
			b = x[bidx]	
			
			#Check Epvs Ord to control
			idx = find_nearest(alphalist, parit["Alpha"])
			if abs(col[aidx]-Epvsordata[aidx, -1]) > col[aidx]:
				continue
			if abs(Epvsordata[aidx, -1] - Epvsordata[aidx, -5]) > 0.001*Epvsordata[aidx, -1]:
				print('Estimator did not converge!')
				continue
			ax[idx].plot(orders, Epvsordata[aidx], label  = r'$\omega_{lowl} = -%g$' %(a))
			ax[idx].plot(orders, Epvsordata[bidx], label  = r'$\omega_{uppl} = -%g$' %(b))
			
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
			
			
								
		if not Eproots:
			os.chdir(os.pardir)
			it+=1
			continue
		
		#check maximum order for estimator convergence
		ordit = 0	
		for Epord in Epvsordata[aidx,:]:
			if abs(Epvsordata[aidx, -1] - Epord) < 0.00001 * Epvsordata[aidx, -1]:
				break
			ordit += 1
		Eporda= orders[ordit]
		ordit = 0	
		for Epord in Epvsordata[bidx,:]:
			if abs(Epvsordata[bidx, -1] - Epord) < 0.00001 * Epvsordata[bidx, -1]:
				break
			ordit += 1
		Epordb= orders[ordit]
		Epord= Eporda if Eporda > Epordb else Epordb
		Epconvord.append([parit["Alpha"], parit["Chemical_Potential"], Epord, orders[-1]])

		#Esti Data TransferEpconvord
		Eprootslist.append([parit["Alpha"], Eproots])
		os.chdir(os.pardir)
		it+=1


#Plot controll Plots		
for pit in range(len(ax)):
	ax[pit].set_yscale('linear')
	ax[pit].set_xlabel('Order')
	ax[pit].set_ylabel(r'$-E_p(<Ord)$')
	ax[pit].legend(loc=4)
	fig[pit].savefig('Epvsord_%.1f.pdf' %(alphalist[pit]))
	fig[pit].clear()
		
#join Eproots		
Eprootslistcomp = [[j, []] for j in alphalist]
for line in Eprootslist:
	idx = find_nearest(alphalist, line[0])
	Eprootslistcomp[idx][1] += line[1]
				
#Mean and Std
lineit = 0	
for line in Eprootslistcomp:
	fig4esti[lineit, 0] = line[0]
	fig4esti[lineit, 1] = np.mean(line[1])
	fig4esti[lineit, 2] = np.std(line[1])/np.sqrt(len(line[1])-1)
	lineit+=1
	
#Data Output
np.savetxt('fig4data', fig4esti)




#join Max convergence order and
#plot some Max order convergence
Epconvordlist = [[j, [], []] for j in alphalist]
for line in Epconvord:
	idx = find_nearest(alphalist, line[0])
	line[1] = fig4esti[idx, 1] + line[1]
	Epconvordlist[idx][1].append(line[1])
	Epconvordlist[idx][2].append(line[2])

for line in Epconvordlist:
	if len(line[1]) > 2 :
		plt.plot(line[1], line[2], 'o')
		plt.title(r'Alpha = %.1f' %(line[0]))
		plt.yscale('linear')
		plt.xlabel(r'$E_p - \mu$')
		plt.ylabel(r'$0.001 Convergence Order$')
		plt.savefig('ConvOrdAlpha%.1f' %(line[0]) + '.pdf')
		plt.clf()
	
Epconvord = np.asarray(Epconvord)
sort = np.argsort(Epconvord[:,0])
Epconvord = Epconvord[sort]
np.savetxt('Epconvord', Epconvord)




#Compare Plot
#plt.errorbar(fig4fit[:, 0], -fig4fit[:, 1], fig4fit[:, 2], fmt='o', label='Fit')
plt.errorbar(fig4esti[:, 0], fig4esti[:, 1], fig4esti[:, 2], fmt='r^', label='Estimator')
#peter = np.loadtxt(peter_loc+ '/fig4_qc'+ str(params['Q_Cutoff']))
#plt.plot(peter[:, 0], peter[:, 1], label="Paper")

def Eprenorm(alpha):
	return alpha/np.sqrt(2)/pi *(1+ 1/params["Impurity_Mass"])*params["Q_Cutoff"]
Eprenplt = np.vectorize(Eprenorm)
xren =np.arange(alphalist[0], alphalist[-1], (alphalist[-1]-alphalist[0])/100)
plt.plot(xren, Eprenplt(xren), label=r'$E_{Ren}$') 

plt.title(r'Q_c = ' + str(params['Q_Cutoff']))
plt.yscale('log')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$-Ep^{MC}$')
plt.legend(loc=4)
plt.savefig('fig4.pdf')
plt.clf()


#Epvlie = interp1d(peter[:,0],peter[:,1])

#plt.plot(fig4fit[:, 0], -(fig4fit[:,1]/Epvlie(fig4fit[:,0])), 'bo', label='Fit')
#plt.plot(fig4esti[:, 0], (fig4esti[:,1]/Epvlie(fig4esti[:,0])), 'r^', label='Estimator')
#plt.axhline(y = np.sqrt(2), color ='g')

#plt.title(r'Q_c = ' + str(params['Q_Cutoff']))
#plt.xlabel(r'$\alpha$')
#plt.ylabel(r'$\frac{E_p^{MC}}{E_p^{Vlietinck}}$')
#plt.legend(loc=4)
#plt.savefig('fig4_comp_Vlie.pdf')
#plt.clf()



#------------------------------------------------------Fig6total

fig6esti = fig4esti
fig6esti[:,1] = Eprenorm(fig4esti[:,0])-fig4esti[:,1]
np.savetxt('fig6data', fig6esti)

Gr = np.genfromtxt(peter_loc + 'Grlarge_Vlie.csv',delimiter=',', skip_header=6)

plt.errorbar(fig6esti[:,0], fig6esti[:,1], fig4esti[:, 2], fmt='r^', label='Estimator')
plt.plot(Gr[:, 0], Gr[:, 1]*np.sqrt(2), label="Grusdt")

plt.title(r'Q_c = ' + str(params['Q_Cutoff']))
plt.yscale('linear')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$Ep^{MC}+ E_{Ren}$')
plt.legend()
plt.savefig('fig6.pdf')
plt.clf()

#------------------------------------------------------Fig6small alpha


Gr = np.genfromtxt(peter_loc + 'Grsmall_Vlie.csv',delimiter=',', skip_header=6)
Grerr = np.genfromtxt(peter_loc +'Grsmall_Vlie_err.csv',delimiter=',', skip_header=6)
Grerr[:,1]-=Gr[:,1]
FeyVar = np.genfromtxt(peter_loc + 'Grsmall_FM.csv',delimiter=',', skip_header=6)
MF = np.genfromtxt(peter_loc + 'Grsmall_MF.csv',delimiter=',', skip_header=6)

plt.errorbar(Gr[:, 0], Gr[:, 1], Grerr[:,1], label="Vlietinck in Grusdt paper")
plt.plot(FeyVar[:, 0], FeyVar[:, 1], label = "FeyVar")
plt.plot(MF[:, 0], MF[:, 1], label = "MF")


plt.errorbar(fig6esti[:,0], fig6esti[:,1], fig4esti[:, 2], fmt='r^', label='Estimator')
plt.yscale('linear')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$E_p + E_{Ren} \, [\frac{c}{\xi}]$')
plt.xlim([0, 1])
plt.ylim([-1,0.7])
plt.legend(loc = 0)
plt.savefig('fig6_small_alpha.pdf')
plt.clf()


