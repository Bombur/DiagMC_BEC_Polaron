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

peter_loc = '/project/theorie/h/H.Guertner/BEC_Polaron/data/without_measreweight'

#sample
data = np.loadtxt('data/all_orders')
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='MEASREWEIGHT')

# Peter
peter = np.loadtxt(peter_loc+'/all_orders')
plt.plot(peter[:, 0], peter[:, 1], label='normal')

plt.xlabel(r'$\tau$')
plt.ylabel(r'$\Sigma (\tau) e^{\mu \tau}$')
plt.xlim([0, 5])
plt.legend()
plt.savefig('vgl_all.pdf')

plt.clf()

#sample
data = np.loadtxt('data/first_order')
plt.errorbar(data[:, 0], data[:, 1], data[:, 2], label='MEASREWEIGTH')

# Peter
peter = np.loadtxt(peter_loc+'/first_order')
plt.plot(peter[:, 0], peter[:, 1], label='normal')

plt.xlabel(r'$\tau$')
plt.ylabel(r'first order')
plt.xlim([0, 5])
plt.legend()
plt.savefig('vgl_first.pdf')

plt.clf()

