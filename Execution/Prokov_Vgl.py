import numpy as np
import matplotlib, math
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

smend =100
peter_loc = '/project/theorie/h/H.Guertner/BEC_Polaron/peter'

#Green
#Prokof'ev
Prok= np.genfromtxt(peter_loc +'/Prok_fig2.csv',delimiter=',', skip_header=6)
plt.plot(Prok[:,0], Prok[:,1], label="Prokof'ev")
#sample
data = np.loadtxt('data/all_orders')
nozero = np.nonzero(data[:,1])
data = data[nozero]
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='Sampled')
plt.yscale('log')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G(0, \tau)$')
plt.legend()
plt.savefig('ana/Prok_Vgl.pdf')
plt.clf()


