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

path= '/home/h/H.Guertner/theorie/BEC_Polaron/'

for folder in os.listdir(os.getcwd()):
	if os.path.isdir(folder) and os.path.exists(folder + "/data/all_orders"):
		os.chdir(folder)
		shutil.copy(path + 'data_ana_fig4.py', os.curdir)
		shutil.copy(path + 'BEC_1st_SEomega.m', os.curdir)
		shutil.copy(path + 'BEC_1st_G0SEG0.m', os.curdir)
		try:
			os.system("python3 data_ana_fig4.py")
		except:
			os.chdir(os.pardir)
			continue
		os.chdir(os.pardir)

