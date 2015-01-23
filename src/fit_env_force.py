#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os.path import expanduser
from scipy.optimize import leastsq
from scipy.ndimage import gaussian_filter

# Models
def lin_mod(par,x):
    k0 = par
    return k0*x

def exp_mod(par,x):
    k1, k2 = par
    return k1*np.exp(-k2*x)

def sig_mod(par,x):
    k3, k4 = par
    return 1/(1+np.exp(k3*x - k4))

# Cost functions
def lin_residuals(par,y,x):
    return y-lin_mod(par,x)

def exp_residuals(par,y,x):
    return y-exp_mod(par,x)

def sig_residuals(par,y,x):
    return y-sig_mod(par,x)

# Import data
figs_path = "../figs/"
file_path = "../data/filtered_SturtPlains_v12.csv"
sp_data = pd.read_csv( expanduser(file_path) )
show_plot = True


# Fix the bugger-up with the tail of the NDVI data
sp_data.loc[2384:2408,'NDVI250X'] = sp_data.loc[2370,'NDVI250X']
# Rough idea to remove trees
yraw = sp_data["NDVI250X"]
ymes = yraw - min(yraw)
if show_plot==True:
    fig = plt.figure()
    plt.plot( yraw, '-', color='lightblue', label='Total')
    plt.plot( ymes, '-', color='blue', lw=2, label='Grass' )
    plt.legend(loc=1)
    plt.xlabel('Days since 1-Jan-2008')
    plt.ylabel('NDVI')
    plt.savefig(figs_path+"ndvi_corrected.pdf")
    plt.close(fig)

# Smooth out soil moisture to get the averaged concurrent point matching the
# inflexion points in the NDVI data
xraw = sp_data["SWC10"]
xmes = gaussian_filter( xraw, 10 )
if show_plot==True:
    fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True)
    ax1.plot( xraw, '-', color='pink' )
    ax1.plot( xmes, linestyle='-', color='red', lw=2)
    ax2.plot( xraw-xmes, '.', color='pink' )
    ax1.set_ylabel(r'$\theta_{s}$', fontsize=18)
    ax2.set_ylabel(r'$\sigma$', fontsize=18)
    ax2.set_xlabel('Days since 1-Jan-2008')
    plt.savefig(figs_path+"swc_smoothed.pdf")
    plt.close(fig)

# From the above two techniques, create a new dataframe for the filtered
# versions of SWC and NDVI
sp_data_new = pd.DataFrame({'SWC10':xmes, 'NDVI250X':ymes})

# Find the inflexion points in the NDVI time-series
yd = [0 for j in range(len(ymes))]
ydd = [0 for j in range(len(ymes))]
for i in range(len(ymes)-1):
    yd[i+1] = ymes[i+1] - ymes[i]
    ydd[i+1] = yd[i+1] - yd[i]
if show_plot==True:
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.plot( ymes, color='purple', lw=2)
    ax2.plot( yd, color='red', lw=2, label="$dx/dt$" )
    ax2.plot( ydd, color='blue', lw=2, label="$d^2x/dt^2$" )
    ax1.set_ylabel('NDVI signal')
    ax2.set_ylabel(r'$x_{t+1}-x_{t}$', fontsize=16)
    ax2.set_xlabel('Days since 1-Jan-2008')
    ax2.axis([1,1.9e3,-6e-3,6e-3])
    ax2.legend(loc=1)
    plt.savefig(figs_path+"swc_turningpoints.pdf")
    plt.close(fig)

# Put these in a dataframe to filter sp_data_new
ydiff = pd.DataFrame({'yd1':yd,'yd2':ydd})
# Tolerance control on what points to extract
tol = 1e-2
upper = max(ydd[20:(len(ydd)-20)])*tol
lower = min(ydd[20:(len(ydd)-20)])*tol
sp_data_filt = sp_data_new.loc[(ydiff['yd1']<upper) & (ydiff['yd1']>lower)]

# Least squares minimization solution to explaining the relationship between
# the extrema of SWC and NDVI
lin0 = [1]
exp0 = [1,1]
sig0 = [1,0]
xobs = sp_data_filt['SWC10']
yobs = sp_data_filt['NDVI250X']
lin_res = leastsq( lin_residuals, lin0, args=(yobs,xobs) )
exp_res = leastsq( exp_residuals, exp0, args=(yobs,xobs) )
sig_res = leastsq( sig_residuals, sig0, args=(yobs,xobs) )

# Save the parameter estimates to a CSV table

# Create vectors for model fits
xs = np.arange(0,0.3,1e-3)
ndvi_lin = lin_mod( lin_res[0], xs )
ndvi_exp = exp_mod( exp_res[0], xs )
ndvi_sig = exp_mod( sig_res[0], xs )

# Plot the results
fig = plt.figure()
plt.plot( sp_data_filt["SWC10"], sp_data_filt["NDVI250X"], 'o', color='black' )
plt.plot( xs, ndvi_lin, linestyle='-', color='red', lw=2, label=r"$k_{0}\theta_{s}$" )
plt.plot( xs, ndvi_exp, linestyle='-', color='blue', lw=2, label=r"$k_{1}\exp(k_{2}\theta_{s})$" )
plt.plot( xs, ndvi_sig, linestyle='-', color='purple', lw=2, label=r"$(1+\exp(k_{3}\theta_{s}-k_{4}))^{-1}$" )
plt.xlabel(r'$\theta_{s 10cm}$', fontsize=18)
plt.ylabel('NDVI')
plt.legend(loc=2)
plt.axis([0,0.25,0,0.5])
plt.savefig(figs_path+"phen_fe_fit.pdf")
