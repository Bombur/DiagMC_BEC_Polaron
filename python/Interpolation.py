#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.interpolate import interp1d

#helper function
def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

def interpoldata(taus, sigma_in, nbins, varname):
	"""
	Interpolation of irreguar spaced data

    Args:
        numpy arrays of taus and SE or G0SE
        int of the new number of bins
        string of the variable name for the plot

    Returns:
        A tuple of taus with equal spacing and the corresponding data
    """
	# Sort out the 0 values
	#nozero = np.nonzero(sigma_in)
	#taus = taus[nozero]
	#sigma_in = sigma_in[nozero]
	
	#Interpolation
	intfunc = interp1d(taus, sigma_in)
	dtau = (taus[-1] - 0)/(nbins+1)
	tau_grid = np.arange(dtau/2, taus[-1], dtau, dtype=float)
	sigma_out = intfunc(tau_grid)
	
	#Plot	
	plt.plot(taus, sigma_in, label='Sampled')
	plt.plot(tau_grid, sigma_out, label = 'Interpolated')
	plt.xlabel(r'$\tau$')
	plt.ylabel(varname)
	plt.legend()
	plt.savefig('ana/Transform/interpol_' + varname + '.pdf')
	plt.clf()

	return tau_grid, sigma_out
    
def find_lastbin(sigma_in, irrel):
	"""
	Sort out irrelevant data

    Args:
        numpy array of SE or G0SE
        int of which exponent of ten is irrelevant

    Returns:
        int of the lastbin
    """
	binstart = round(sigma_in.size*0.75)
	lastbin= binstart + find_nearest(sigma_in[binstart:], 10**(irrel))
	lastbin2= binstart + find_nearest(sigma_in[binstart:], 0)
	if lastbin2 < lastbin:
		lastbin=lastbin2
		
	if lastbin > sigma_in.size -5 :
		return sigma_in.size -5
	elif abs(sigma_in[lastbin]- 10**(irrel)) > 1:
		return sigma_in.size-10
	else:
		return lastbin	
 
if __name__ == "__main__":
	data = np.loadtxt('../data/zero_order')
	lastbin = find_lastbin(data[:,1], -20)
	print(lastbin)
	_, dataout= interpoldata(data[:lastbin,0], data[:lastbin,1], 200, 'G0SE')

    
