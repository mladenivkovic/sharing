#!/usr/bin/python3

#====================================================================
# Plots the results of eval_galaxies.f90 scripts.
# This script is called by the eval_galaxies.sh script, but you can 
# also call it manually:
#   $ plot_fortran_correlation.py output_XXXXX 
# it looks for the following files:
# Pk_all.txt, Pk_sub.txt, correlation_all.txt, correlation_sub.txt
#====================================================================

import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os.path import exists
from time import sleep




#------------------------
# Get case from cmdline
#------------------------


srcdir = argv[1]



pkfile_all = srcdir+'/Pk_all.txt'
pkfile_sub = srcdir+'/Pk_sub.txt'
xifile_all = srcdir+'/correlation_all.txt'
xifile_sub = srcdir+'/correlation_sub.txt'

if not exists(srcdir):
    print("Didn't find directory ", srcdir, "that you passed as argument. Aborting")
    quit()

filelist=[pkfile_all, pkfile_sub, xifile_all, xifile_sub]
for f in filelist:
    if not exists(f):
        print("Didn't find file ", f)
        quit(2)


h = 0.704


def obs_correlation(r):
    #  https://arxiv.org/pdf/astro-ph/0301280.pdf
    return (r/(5.77/h))**(-1.80)

#============================
def read_power_spectrum():
#============================
    """
    Reads in power spectrum from observations from file.
    """
    k_eff, k_low, k_high, P, dP = np.loadtxt('/home/mivkov/UZH/Masterarbeit/masterarbeit/observational_data/2dF_power_spectrum.txt',
        skiprows=4, unpack=True)

    k_high -= k_eff
    k_low = k_eff-k_low


    # change units from h Mpc-1 to Mpc-1
    
    k_high  /= h
    k_low   /= h
    k_eff   /= h
    P       *= h**3
    dP      *= h**3

    return k_eff, (k_low, k_high), P, dP




#-------------------------------
# Read and plot data
#-------------------------------

k_sub, Pk_sub = np.loadtxt(pkfile_sub, dtype='float', unpack=True, skiprows=2)
k_all, Pk_all = np.loadtxt(pkfile_all, dtype='float', unpack=True, skiprows=2)
r_sub, xi_sub = np.loadtxt(xifile_sub, dtype='float', unpack=True, skiprows=2)
r_all, xi_all = np.loadtxt(xifile_all, dtype='float', unpack=True, skiprows=2)

fig = plt.figure(figsize=(16,10))

ax1 = fig.add_subplot(121)
ax1.loglog(k_sub, Pk_sub, label='without orphans')
ax1.loglog(k_all, Pk_all, label='with orphans')

k_obs, k_err, P_obs, P_err = read_power_spectrum()
ax1.errorbar(k_obs, P_obs, yerr=P_err, xerr=k_err, label='observational data')
ax1.set_title("Power Spectrum")
ax1.set_xlabel(r"$k$")
ax1.set_ylabel(r"$P(k)$")
ax1.legend()


ax2 = fig.add_subplot(122)
ax2.semilogy(r_sub, xi_sub, label='without orphans')
ax2.semilogy(r_all, xi_all, label='with orphans')
ax2.semilogy(r_all, obs_correlation(r_all), label='observational data')

ax2.set_title("Correlation")
ax2.set_xlabel(r"$r$")
ax2.set_ylabel(r"$\xi(r)$")
ax2.legend()

plt.savefig(srcdir+"/correlation.png")

