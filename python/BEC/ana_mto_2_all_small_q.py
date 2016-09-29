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

taus = [i for i in params['Taus']]
taumax = taus[-1]

# data for the comparison    
peter_loc = params["Peter_Path"]
scale= 'log'
    
# cleanup
def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured


# definitions
pi = math.pi


#--------------- plot Q Statistics
#sample
for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/data/all_orders"):
		os.chdir(folder)
		data = np.loadtxt('data/stats/qs_tot', skiprows=1)
		norm = data[:,1].sum()
		plt.plot(data[:, 0], data[:, 1]/norm, label = folder)
		os.chdir(os.pardir)
plt.xlabel(r'Phonon Momentum $q$')
plt.ylabel(r'$P(q)$ [%]')
plt.xlim([0, 200])
plt.legend()
plt.savefig('plot_q.pdf')

plt.clf()

#--------- plot Green-function all

for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/data/all_orders"):
		os.chdir(folder)
		minmax= np.loadtxt('data/stats/minmax', skiprows = 3)
		tmp = "%i" %(minmax[-1,1])
		data = np.loadtxt('data/all_orders')
		plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label= folder + " MO_" +tmp)
		os.chdir(os.pardir)

# Peter
peter = np.loadtxt(peter_loc+'/fig3total')
plt.plot(peter[:, 0], peter[:, 1], label= "fig3total")

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma (\tau)$')
plt.xlim([0, taumax])
#plt.ylim([1e-5, 1e5])
#plt.legend()
plt.savefig('plot_all.pdf')

plt.clf()

#--------- plot Green-function all small

for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/data/all_orders"):
		os.chdir(folder)
		minmax= np.loadtxt('data/stats/minmax', skiprows = 3)
		tmp = "%i" %(minmax[-1,1])
		data = np.loadtxt('data/all_orders')
		plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label= folder + " MO_" +tmp)
		os.chdir(os.pardir)

# Peter
peter = np.loadtxt(peter_loc+'/fig3small')
plt.plot(peter[:, 0], peter[:, 1], label= "fig3small")

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma (\tau)$')
plt.xlim([0, 1e-3])
plt.ylim([1e2, 1e6])
plt.legend(loc=4)
plt.savefig('plot_all_small.pdf')

plt.clf()

#--------- plot 2nd order

for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/data/all_orders"):
		os.chdir(folder)
		data = np.loadtxt('data/second_ordera')
		plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label= folder)
		os.chdir(os.pardir)

# Peter
peter = np.loadtxt(peter_loc+'/fig3small')
plt.plot(peter[:, 0], peter[:, 1], label= "fig3small")

plt.yscale(scale)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma (\tau)$')
plt.xlim([0, 1e-3])
plt.ylim([1e2, 1e5])
plt.legend(loc=4)
plt.savefig('plot_2_small.pdf')

plt.clf()




