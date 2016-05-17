import numpy as np
import matplotlib, math
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

smend =100
peter_loc = '/project/theorie/h/H.Guertner/BEC_Polaron/peter'

def Eprenorm(alpha, qc):
	return alpha/np.sqrt(2)/math.pi *(1+ 1/0.263158)*qc

#-------------------------------SE >1
#Vlietinck
vliesm = np.loadtxt(peter_loc + '/fig3small')
vliesm2up= np.genfromtxt(peter_loc +'/fig3small2up.csv',delimiter=',', skip_header=6)
vliesm2down= np.genfromtxt(peter_loc +'/fig3small2down.csv',delimiter=',', skip_header=6)
vlie = np.loadtxt(peter_loc + '/fig3big')
print(-790 , '\t', Eprenorm(5,200)/np.sqrt(2))
plt.plot(vlie[:,0], vlie[:,1]*np.sqrt(2)*np.exp(-(-790 + (Eprenorm(5,200)/np.sqrt(2)))*vlie[:,0]), label="Vlietinck")
#sample
data = np.loadtxt('data/SE/SE_>1')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1]*np.exp(-(-55 + Eprenorm(5,10))*data[:,0]), data[:, 2], label='Sampled')
plt.yscale('log')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma^{>1}(0, \tau)$')
plt.legend()
plt.savefig('ana/SE/SE_>1_Vlie.pdf')
plt.clf()

#small plot
plt.plot(vliesm[:,0], vliesm[:,1]*np.sqrt(2), label="Vlietinck 1")
plt.plot(vliesm2up[:,0], vliesm2up[:,1]*np.sqrt(2)/10, label="Vlietinck Upper Limit")
plt.plot(vliesm2down[:,0], vliesm2down[:,1]*np.sqrt(2)/10, label="Vlietinck Lower Limit")
plt.errorbar(data[:smend,0],data[:smend,1], data[:smend,2], label="Sampled")
plt.yscale('log')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma^{>1}(0, \tau)$')
plt.legend()
plt.savefig('ana/SE/SE_>1_Vlie_small.pdf')
plt.clf()

