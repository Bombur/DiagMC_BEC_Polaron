#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.special import dawsn
import Interpolation as Ip
import re

#define constants
pi = np.pi

# import parameters
with open('DiagMC_BEC.json') as params_file:
    params = json.load(params_file)

mu = params['Chemical_Potential']
p_initial = params['Momentum']
relm= params['Impurity_Mass']
bins = params['Bins']
tau_equi = False if len(bins)>2 else True	
irrel = params['Irrelevance_Order']
intbin = params['Interpol_Bin']

#Read in complex data
def cfcomp(x):
    m = re.findall('[-\d.]+', x.decode("utf-8"))
    return np.complex(float(m[0]),float(m[1]))
    
def mat_cfcomp(x):
    m = re.findall('[-\d.]+', x.decode("utf-8"))
    return np.complex(float(m[0]),float(m[1]))

def cfreal(x):
    m = re.findall('[-\d.]+', x.decode("utf-8"))
    return float(m[0])
    

def green_tau_vac(tau, p):
    """
    Vacuum Green's function in imaginary-time-momentum representation
    """
    if params['BEC_Polaron']:
        return np.exp(-(p**2 * 0.5/relm * np.sqrt(2) - mu)*tau)
    else:
	    return np.exp(-(p**2 * 0.5 - mu)*tau)

def green_omega_vac(omega, p):
    """
    Vacuum Green's function in frequency-momentum representation
    """
    if params['BEC_Polaron']:
        return -1./(1j*omega - (p**2 * 0.5/relm * np.sqrt(2)) + mu)
    else:
        return -1./(1j*omega - p**2 * 0.5 + mu)


def fourier_to_omega_fft(sigma_in, dtau):
    """
    Fast-Fourier-transformation of measured self-energy

    Args:
        sigma_in: NumPy-array to be Fourier transformed. It is assumed real,
                    with data at dtau/2, 3*dtau/2, ...
        dtau: distance in tau between two values of sigma_in

    Returns:
        A tuple of omega-data and the corresponding Fourier transformed values
    """
    sigma_out = np.fft.rfft(sigma_in)
    # Fourier convention of NumPy is different from many-body:
    sigma_out = dtau * sigma_out.conjugate()

    # get corresponding omega-values
    omega_out = np.arange(sigma_out.size, dtype=float)
    omega_out *= 2.*pi / dtau / sigma_in.size

    # consider offset
    sigma_out *= np.exp(1j * dtau/2. * omega_out)
    return omega_out, sigma_out


def fourier_to_omega(sigma_in, dtau):
    """
    Fourier-transformation of measured self-energy

    Args:
        sigma_in: NumPy-array to be Fourier transformed. It is assumed real,
                    with data at dtau/2, 3*dtau/2, ...
        dtau: distance in tau between two values of sigma_in

    Returns:
        A tuple of omega-data and the corresponding Fourier transformed values
    """
    sigma_out = np.zeros(sigma_in.size, dtype=complex)
    omega_out = np.arange(sigma_out.size, dtype=float)
    omega_out *= pi/ dtau / sigma_out.size
    for i, omega in enumerate(omega_out):
        sigma_out[i] = 0
        for j, sigma in enumerate(sigma_in):
            current_tau = dtau/2. + j * dtau
            sigma_out[i] += sigma * np.exp(1j * omega * current_tau)
    sigma_out *= dtau
    return omega_out, sigma_out


def fourier_to_tau_fft(green_in, domega):
    """
    Fast-Fourier-transformation to imaginary time Green's function

    Args:
        green_in: NumPy-array to be Fourier transformed. It is assumed hermitian
                    symmetric, with data at 0, domega, 2*domega ...
        domega: distance in omega between two values of green_in

    Returns:
        A tuple of omega-data and the corresponding Fourier transformed values
    """
    green_out = np.fft.irfft(green_in.conjugate())
    green_out *= domega * green_out.size/2./pi

    tau_out = np.arange(green_out.size, dtype=float)
    tau_out *= 2.*pi / domega / green_out.size
    return tau_out, green_out


def fourier_to_tau(green_in, domega):
    """
    Fourier-transformation to imaginary time Green's function

    Args:
        green_in: NumPy-array to be Fourier transformed. It is assumed hermitian
                    symmetric, with data at 0, domega, 2*domega ...
        domega: distance in omega between two values of green_in

    Returns:
        A tuple of omega-data and the corresponding Fourier transformed values
    """
    green_out = np.zeros(green_in.size, dtype=float)
    tau_out = np.arange(green_out.size, dtype=float)
    tau_out *= pi / domega / green_out.size
    for i, tau in enumerate(tau_out):
        green_out[i] = 0
        for j, green in enumerate(green_in):
            current_omega = j * domega
            if j == 0:
                green_out[i] += (green * np.exp(-1j * tau * current_omega)).real
            else:
                green_out[i] += (2. * green * np.exp(-1j * tau * current_omega)).real
    green_out *= domega/2./pi
    return tau_out, green_out


def dyson_g0se(green_omega_vac, sigma_in):
    """
    Apply Dyson's equation to input vacuum Green's function
    and self-energy (multiplied with vacuum Green's function)

    Args:
        green_omega_vac: vacuum Green's function data at desired omega-points
        sigma_in: Sigma-G0 data at desired omega-points

    Returns:
        tau grid and Interacting Green's function
    """
    return ((1. - sigma_in)/green_omega_vac)**-1
    
    
def dyson(green_omega_vac, sigma_in):
    """

    Apply Dyson's equation to input vacuum Green's function
    and self-energy (multiplied with vacuum Green's function)

    Args:
        green_omega_vac: vacuum Green's function data at desired omega-points
        sigma_in: Sigma data at desired omega-points


    Returns:
        tau grid and Interacting Green's function
    """
    return (1./green_omega_vac - sigma_in)**-1


def green_trafo():
    """
    Calculate and plot the imaginary-time Green's function out of selfenergy, G0SE and G0SE in Matsubara
    """
    tau_grid = np.loadtxt('data/all_orders')[:, 0]

    #Sigma Transformation
    sigma_tau = np.loadtxt('data/SE/SE_>1')[:, 1]
    sigma_tau_first = np.loadtxt('data/SE/mat_1st_se')[:, 1]
    if not tau_equi:
        lastbin = Ip.find_lastbin(sigma_tau, irrel)
        tau_grid_SE, sigma_tau = Ip.interpoldata(tau_grid[:lastbin], sigma_tau[:lastbin], intbin, 'SE')
        _, sigma_tau_first = Ip.interpoldata(tau_grid[:lastbin], sigma_tau_first[:lastbin], intbin, 'SE1')
    else:
    	tau_grid_SE = tau_grid
    omega_grid_SE, sigma_omega = fourier_to_omega_fft(sigma_tau, tau_grid_SE[2] - tau_grid_SE[1])
    _, sigma_omega_first = fourier_to_omega_fft(sigma_tau_first, tau_grid_SE[2] - tau_grid_SE[1])
    sigma_omega += sigma_omega_first


    if params['Self_Energy']:  
        #G0Sigma Transformation
        sigma_G0_tau = np.loadtxt('data/all_orders')[:, 1]
        sigma_G0_tau_first = np.loadtxt('data/mat_1st_g0se')[:, 1]
        if not tau_equi:
            lastbin = Ip.find_lastbin(sigma_G0_tau, irrel)
            tau_grid_G0SE, sigma_G0_tau = Ip.interpoldata(tau_grid[:lastbin], sigma_G0_tau[:lastbin], intbin, 'G0SE')
            _, sigma_G0_tau_first = Ip.interpoldata(tau_grid[:lastbin], sigma_G0_tau_first[:lastbin], intbin, 'G0SE1')
        else:
    	    tau_grid_G0SE = tau_grid
    	#FFT
        omega_grid_G0SE, sigma_G0_omega = fourier_to_omega_fft(sigma_G0_tau, tau_grid_G0SE[2] - tau_grid_G0SE[1])
        _, sigma_G0_omega_first = fourier_to_omega_fft(sigma_G0_tau_first, tau_grid_G0SE[2] - tau_grid_G0SE[1])
        sigma_G0_omega += sigma_G0_omega_first
        plt.plot(omega_grid_G0SE, sigma_G0_omega.real, label='FFT Real')
        plt.plot(omega_grid_G0SE, sigma_G0_omega.imag, label='FFT Imag')
		
		##Brute Force
        omega_grid_G0SE2, sigma_G0_omega2 = fourier_to_omega(sigma_G0_tau, tau_grid_G0SE[2] - tau_grid_G0SE[1])
        _, sigma_G0_omega_first2 = fourier_to_omega(sigma_G0_tau_first, tau_grid_G0SE[2] - tau_grid_G0SE[1])
        sigma_G0_omega2 += sigma_G0_omega_first2
        plt.plot(omega_grid_G0SE2, sigma_G0_omega2.real, label='Brute Real')
        plt.plot(omega_grid_G0SE2, sigma_G0_omega2.imag, label='Brute Imag')
    
    
    
    #Sampled G0SEiw
    G0SEiwdata = np.genfromtxt('data/G0SEiw/G0SEiw_all', converters = {0: cfreal, 1: cfcomp, 2:cfcomp})
    omega_grid_G0SEiw = G0SEiwdata['f0']
    G0SEiw = G0SEiwdata['f1']/abs(mu)
    
    #G0SEiw_first = np.genfromtxt('data/G0SEiw/mat_1st_g0seiw', converters = {0:cfcomp})
    #plt.plot(omega_grid_G0SEiw, G0SEiw_first.real, label='Mat Real part')
    #plt.plot(omega_grid_G0SEiw, G0SEiw_first.imag, label='Mat Imaginary part')
    
    
    plt.plot(omega_grid_G0SEiw, G0SEiw.real, label='Samp Real part')
    plt.plot(omega_grid_G0SEiw, G0SEiw.imag, label='Samp Imaginary part')
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$\Sigma G_{0}$')
    plt.legend()
    plt.savefig('ana/Transform/g0selfenergy_omega.pdf')
    plt.clf()

    plt.plot(omega_grid_SE, sigma_omega.real, label='Real part')
    plt.plot(omega_grid_SE, sigma_omega.imag, label='Imaginary part')
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$\Sigma$')
    plt.legend()
    plt.savefig('ana/Transform/selfenergy_omega.pdf')
    plt.clf()


    green_omega_vac_SE = green_omega_vac(omega_grid_SE, p_initial)
    green_omega_SE = dyson(green_omega_vac_SE, sigma_omega)
    
    green_omega_vac_G0SEiw = green_omega_vac(omega_grid_G0SEiw, p_initial)
    green_omega_G0SEiw = dyson_g0se(green_omega_vac_G0SEiw, G0SEiw)
    
    if params['Self_Energy']:
        green_omega_vac_G0SE = green_omega_vac(omega_grid_G0SE, p_initial)
        green_omega_G0SE = dyson_g0se(green_omega_vac_G0SE, sigma_G0_omega)

        plt.plot(omega_grid_G0SE, green_omega_G0SE.real, label=r'$G_0\Sigma$ Re')
        plt.plot(omega_grid_G0SE, green_omega_G0SE.imag, label=r'$G_0\Sigma$ Im')
    plt.plot(omega_grid_SE, green_omega_SE.real, label=r'$\Sigma$ Re')
    plt.plot(omega_grid_SE, green_omega_SE.imag, label=r'$\Sigma$ Im')
    plt.plot(omega_grid_G0SEiw, green_omega_G0SEiw.real, label=r'Samp $G_0\Sigma$ Re')
    plt.plot(omega_grid_G0SEiw, green_omega_G0SEiw.imag, label=r'Samp $G_0\Sigma$ Im')
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'G')
    plt.legend()
    plt.savefig('ana/Transform/green_omega.pdf')
    plt.clf()

    # subtract high-frequency tail
    green_omega_SE -= green_omega_vac_SE
    green_omega_G0SEiw -= green_omega_vac_G0SEiw
    if params['Self_Energy']:
        green_omega_G0SE -= green_omega_vac_G0SE

        plt.plot(omega_grid_G0SE, green_omega_G0SE.real, label=r'$G_0\Sigma$ Re')
        plt.plot(omega_grid_G0SE, green_omega_G0SE.imag, label=r'$G_0\Sigma$ Im')
    plt.plot(omega_grid_SE, green_omega_SE.real, label=r'$\Sigma$ Re')
    plt.plot(omega_grid_SE, green_omega_SE.imag, label=r'$\Sigma$ Im')
    plt.plot(omega_grid_G0SEiw, green_omega_G0SEiw.real, label=r'Samp $G_0\Sigma$ Re')
    plt.plot(omega_grid_G0SEiw, green_omega_G0SEiw.imag, label=r'Samp $G_0\Sigma$ Im')
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$G-G_0$')
    plt.legend()
    plt.savefig('ana/Transform/green_omega_reduced.pdf')
    plt.clf()

    plt.plot(omega_grid_SE[1:], np.log(abs(green_omega_SE[1:]))/np.log(omega_grid_SE[1:]), label=r'$\Sigma$' ) 
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$\log(|G-G_0|)$')
    plt.legend()
    plt.savefig('ana/Transform/green_omega_reduced_log.pdf')
    plt.clf()

    tau_grid_SE, green_tau_SE = fourier_to_tau_fft(green_omega_SE, omega_grid_SE[1])
    tau_grid_G0SEiw, green_tau_G0SEiw = fourier_to_tau_fft(green_omega_G0SEiw, omega_grid_G0SEiw[1])
    if params['Self_Energy']:
        tau_grid_G0SE, green_tau_G0SE = fourier_to_tau_fft(green_omega_G0SE, omega_grid_G0SE[1])

    # add high-frequency tail
        green_tau_G0SE += green_tau_vac(tau_grid_G0SE, p_initial)
        plt.plot(tau_grid_G0SE, green_tau_G0SE, label=r'$G_0\Sigma$')
        
    green_tau_SE += green_tau_vac(tau_grid_SE, p_initial)
    green_tau_G0SEiw += green_tau_vac(tau_grid_G0SEiw, p_initial)

    
    plt.plot(tau_grid_SE, green_tau_SE, label=r'$\Sigma$')
    plt.plot(tau_grid_G0SEiw, green_tau_G0SEiw, label=r'Samp $G_0\Sigma$')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$G$')
    plt.legend()
    plt.savefig('ana/Transform/green_tau.pdf')
    plt.clf()
	
    if params['Self_Energy']:
        return tau_grid_SE, green_tau_SE, tau_grid_G0SEiw, green_tau_G0SEiw, tau_grid_G0SE, green_tau_G0SE
    else:
        return tau_grid_SE, green_tau_SE, tau_grid_G0SEiw, green_tau_G0SEiw


if __name__ == "__main__":
    # run tests

    # TEST 1: Tau -> Omega
    # --------------------

    # define sample data to be transformed
    def test_function_tau(tau):
        return np.exp(-tau**2)

    def analytic_result_omega(omega):
        return 0.5*np.exp(-omega**2/4.) * np.sqrt(pi) + 1j * dawsn(omega/2)

    tau = np.linspace(0, 10, 100) + 0.05
    data_tau = test_function_tau(tau)

    omega_out1, sigma_out1 = fourier_to_omega_fft(data_tau, 0.1)
    omega_out2, sigma_out2 = fourier_to_omega(data_tau, 0.1)

    # plot results
    fig = plt.figure()
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax2 = fig.add_axes([0.3, 0.33, 0.58, 0.35])
    ax1.plot(omega_out1, sigma_out1.real, label='FFT real')
    ax1.plot(omega_out1, sigma_out1.imag, label='FFT imaginary')
    ax1.plot(omega_out2, sigma_out2.real, label='Brute force real')
    ax1.plot(omega_out2, sigma_out2.imag, label='Brute force imaginary')
    omega_range = np.linspace(0, 50, 1000)
    ax1.plot(omega_range, analytic_result_omega(omega_range).real, label='Analytic real')
    ax1.plot(omega_range, analytic_result_omega(omega_range).imag, label='Analytic imaginary')
    ax1.set_xlabel(r'$\omega$')
    ax1.set_xlim(0, 40)
    ax1.legend(ncol=2)

    # inset
    func_rel1 = sigma_out1 - analytic_result_omega(omega_out1)
    func_rel2 = sigma_out2 - analytic_result_omega(omega_out2)
    ax2.plot(omega_out1, func_rel1.real, label=r'$\Delta$(FFT - analytic)(real)')
    ax2.plot(omega_out1, func_rel1.imag, label=r'$\Delta$(FFT - analytic)(imaginary)')
    ax2.plot(omega_out2, func_rel2.real, label=r'$\Delta$(Brute force - analytic)(real)')
    ax2.plot(omega_out2, func_rel2.imag, label=r'$\Delta$(Brute force - analytic)(real)')
    ax2.set_xlabel(r'$\omega$')
    ax2.legend(prop={'size':6})

    plt.savefig('ana/Transform/plot_fourier_omega_test.pdf')
    plt.clf()

    # TEST 2: Omega -> Tau
    # --------------------

    # define sample data to be transformed
    def test_function_omega(omega):
        return 1./(omega - 2j)/(omega + 1j)

    def analytic_result_tau(tau):
        return 1./3.*np.exp(-tau)

    omega_range = np.linspace(0, 30, 100)
    data_omega = test_function_omega(omega_range)
    tau_out1, green_out1 = fourier_to_tau_fft(data_omega, omega_range[1])
    tau_out2, green_out2 = fourier_to_tau(data_omega, omega_range[1])

    # plot results
    fig = plt.figure()
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax2 = fig.add_axes([0.3, 0.35, 0.58, 0.35])
    ax1.plot(tau_out1, green_out1, label='FFT')
    ax1.plot(tau_out2, green_out2, label='Brute Force')

    ax1.plot(tau, analytic_result_tau(tau), label='Analytic')
    ax1.set_xlim(0, 10)
    ax1.set_xlabel(r'$\tau$')
    ax1.legend()

    # inset
    ax2.plot(tau_out1, green_out1 - analytic_result_tau(tau_out1), label=r'$\Delta$(FFT - analytic)')
    ax2.plot(tau_out2, green_out2 - analytic_result_tau(tau_out2), label=r'$\Delta$(Brute force - analytic)')
    ax2.set_xlabel(r'$\tau$')
    ax2.legend()

    plt.savefig('ana/Transform/plot_fourier_tau_test.pdf')
    plt.clf()
    
    green_trafo()




