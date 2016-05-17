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
		shutil.copy(path + 'fig4.py', os.curdir)

		try:
			os.system("python3 fig4.py")
		except:
			os.chdir(os.pardir)
			continue
		os.chdir(os.pardir)

