#!/usr/bin/env python3
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, errno, sys
import math
import json
import shutil
#from scipy.optimize import curve_fit
from scipy import optimize
from scipy.interpolate import interp1d, splev, splrep

path= '/home/h/H.Guertner/theorie/BEC_Polaron/python/BEC/'

for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder): 
		os.chdir(folder)
		for subfolder in os.listdir(os.getcwd()):
			if os.path.isdir(subfolder) and os.path.exists(subfolder + "/data/all_orders"):
				os.chdir(subfolder)
				allord = np.loadtxt('data/all_orders')
				print(allord)
				allord[:,1:] *= 0.05
				#print(allord)
				#np.savetxt('data/all_orders2', allord)
				Epvsws = np.loadtxt('data/Ep/Epvsws')
				#print(Epvsws)
				Epvsws[:,1:] *= 0.05
				#print(Epvsws)
				#np.savetxt('data/Ep/Epvsws2', Epvsws)
				Epvsord = np.loadtxt('data/Ep/Epvsord')
				#print(Epvsord)
				Epvsord[:,1:] *= 0.05
				#print(Epvsord)
				#np.savetxt('data/Ep/Epvsord2', Epvsord)
				
				#try:
					#os.system("python3 data_ana_fig4.py")
				#except:
				#	os.chdir(os.pardir)
				#	continue
				os.chdir(os.pardir)

		os.chdir(os.pardir)

