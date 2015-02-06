#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from os import listdir
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter
from numpy import linalg as la

# Functions
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

def min_chi2(fun,y,x):
    return lambda par: ((fun(par,x)-y)**2).sum()

def get_errors(opt_res, Nx):
    hessian = opt_res['hess_inv']
    fmin0 = opt_res['fun']
    par = opt_res['x']
    Mp = len(par)
    H_sol = np.diag(la.solve(hessian,np.eye(Mp)))
    coeff = Mp*fmin0/(Nx-Mp)
    return np.sqrt(coeff*H_sol)

# Plot functions
def plot_smooth_swc(xraw, xmes):
    fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True)
    ax1.plot( xraw, '-', color='pink' )
    ax1.plot( xmes, linestyle='-', color='red', lw=2)
    ax2.plot( xraw-xmes, '.', color='pink' )
    ax1.set_ylabel(r'$\theta_{s}$', fontsize=18)
    ax2.set_ylabel(r'$\sigma$', fontsize=18)
    ax2.set_xlabel('Days since 1-Jan-2008')
#    plt.savefig(opath+"swc_smoothed.pdf")
#    plt.close(fig)

def plot_grass_predict(df, yraw, ymes):
    yraw = df["NDVI250X"]
    ymes = df["NDVI250X_Grass"]
    plt.plot( yraw, '-', color='lightblue', label='Total')
    plt.plot( ymes, '-', color='blue', lw=2, label='Grass' )
    plt.legend(loc=1)
    plt.xlabel('Days since 1-Jan-2008')
    plt.ylabel('NDVI')
#    plt.savefig(opath+"ndvi_corrected.pdf")
#    plt.close(fig)

def plot_inflexion_points(ymes, yd, ydd):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.plot( ymes, color='purple', lw=2)
    ax2.plot( yd, color='red', lw=2, label="$dx/dt$" )
    ax2.plot( ydd, color='blue', lw=2, label="$d^2x/dt^2$" )
    ax1.set_ylabel('NDVI signal')
    ax2.set_ylabel(r'$x_{t+1}-x_{t}$', fontsize=16)
    ax2.set_xlabel('Time-series (t)')
    ax2.axis([1,len(ymes),-6e-3,6e-3])
    ax2.legend(loc=1)
    plt.show()
    #plt.savefig(opath+"swc_turningpoints.pdf")
    #plt.close(fig)

def plot_force_optfits(dataset, p_table, \
                       xlabel="SWC_smooth", ylabel="NDVI_grass"):
    # Create vectors for model fits
    xs = np.arange(0,0.3,1e-3)
    ndvi_lin = lin_mod( p_table[p_table["lin"]]['value'], xs )
    ndvi_exp = exp_mod( p_table[p_table["exp"]]['value'], xs )
    ndvi_sig = sig_mod( p_table[p_table["sig"]]['value'], xs )

    # Plot the results
    plt.plot( dataset[xlabel], dataset[ylabel], 'o', color='black' )
    plt.plot( xs, ndvi_lin, linestyle='-', color='red', lw=2, label=r"$k_{0}\theta_{s}$" )
    plt.plot( xs, ndvi_exp, linestyle='-', color='blue', lw=2, label=r"$k_{1}\exp(k_{2}\theta_{s})$" )
    plt.plot( xs, ndvi_sig, linestyle='-', color='purple', lw=2, label=r"$(1+\exp(k_{3}\theta_{s}-k_{4}))^{-1}$" )
    plt.xlabel(r'$\theta_{s 10cm}$', fontsize=18)
    plt.ylabel('NDVI')
    plt.legend(loc=2)
    #plt.axis([0,0.25,0,0.5]) #plt.savefig(opath+"_phen_fe_fit.pdf")
    #plt.close()
    plt.show()

def build_report():
    return None

def get_all_site_names(ec_files):
    file_name = [ re.compile('^\w+').findall(f) for f in ec_files ]
    split_name = [ re.split('_',f) for f in sum(file_name,[]) ]
    name_ver = [ sublist[-2] for sublist in split_name ]
    return name_ver

# Main
def import_data(dat_path, includes=r'^filtered.*'):
    """
    Function reads in datasets stored in the /data folder in the upper
    directory, attaches site names based on a regex search on the file names
    and concatenates the list of datasets into one Pandas dataframe for
    analysis
    """
    # get file names
    file_names = [ f for f in listdir(dat_path) if re.match(includes,f) ]
    # read files in
    datasets = [ pd.read_csv(dat_path+x) for x in file_names ]
    # find site names
    site_names = get_all_site_names( file_names )
    # attach site names relating to each imported dataset
    for i, df in enumerate(datasets):
        df['Site'] = site_names[i]
    # concatenate and return final dataframe
    return datasets

def grass_correct_data(dataset, xvar="SWC10", yvar="NDVI250X"):
    """
    Function filters data to enable a better approximation of the
    environmental forcing that drives the momentum of the pendulum. In this
    case we have smoothed the soil moisture response to reduce the number of
    extrema points and better align with the extrema in the NDVI time-series
    """
    x_raw = dataset[xvar]
    y_raw = dataset[yvar]
    # smooth SMC to obtain more identifiable extrema points (removes noise)
    x_approx = gaussian_filter(x_raw, 10)
    # adjust NDVI by subtracting minimum time-series value to approximate grasses
    y_approx = y_raw - min(y_raw)
    # assign new values as extra columns to the dataframe
    dataset['SWC_smooth'] = x_approx
    dataset['NDVI_grass'] = y_approx
    # return to user
    return dataset

def find_ts_extrema(dataset, var="NDVI_grass"):
    """
    Finds the maxima and minima or turning points of a time-series
    """
    # Find the inflexion points in the NDVI time-series
    dataset['dy1/dt1'] = dataset[var].shift(1).fillna(0)-dataset[var]
    dataset['dy2/dt2'] = dataset['dy1/dt1'].shift(1).fillna(0)-dataset['dy1/dt1']
    return dataset

def get_extrema_points(dataset, tol=1e-2):
    """
    Extract data points along the dy/dx time-series that intersects with the
    dy2/dx2 time-series. This should approximate points that occur around the
    extrema of the dependent variable time-series
    """
    # Put these in a dataframe to filter sp_data_new
    upper = max(dataset["dy2/dt2"].ix[2:])*tol
    lower = min(dataset["dy2/dt2"].ix[2:])*tol
    dataset_ex = dataset.loc[(dataset['dy1/dt1']<upper) & (dataset['dy1/dt1']>lower)]
    if len(dataset_ex.index) <= 1:
        print "ERROR: Tolerance is too high :: not enough data to optimize on"
        return None
    else:
        return dataset_ex

def df_pop_site(dseries):
    """
    Return unique key values from dictionary
    """
    return [ s for s in set(dseries) ]

def create_label(dseries):
    """
    If there is more than one unique label value, then take the first two
    characters and join them into a new string with an underscore linker.
    Otherwise just return the singular name as the label.
    """
    sites = df_pop_site(dseries)
    if len(sites) > 1:
        prefix = [ char[:2] for char in sites ]
        return '_'.join(prefix)
    else:
        return sites.pop()

def approx_optim_force(extrema_df, f_mod, p0=[1,2], xlabel="SWC_smooth", ylabel="NDVI_grass"):
    """
    This function acts as wrapper to fit some arbitrary univariate model given
    X and Y data series. Some dataframe is passed and the two variables of
    interest are extracted based on the two label values, and then optimised
    on. Returns tabular dataframe giving parameter estimates and their errors.
    """
    xobs = extrema_df[xlabel]
    yobs = extrema_df[ylabel]
    sig_res = minimize( min_chi2(f_mod, yobs, xobs), p0 )
    sig_err = get_errors(sig_res, len(yobs))
    site_lab = create_label(extrema_df["Site"])
    return pd.DataFrame({'site':site_lab, 'value':sig_res['x'], 'error':sig_err})

def optimize_all_comb(dataset):
    # ensemble optimisation
    all_data = pd.concat(dataset)
    all_res = approx_optim_force( all_data, sig_mod )
    # out-of-sample optimisation
    sites = df_pop_site(all_data["Site"])
    out_res = []
    for i in sites:
        out_res.append( approx_optim_force( all_data.query("Site != 'i'"), sig_mod ) )
    # in-sample optimisation
    ind_res = pd.concat( map( lambda x: approx_optim_force(x, sig_mod), dataset ) )
    # join all tables of optimisation combinations and return
    return pd.concat([all_res,ind_res]+out_res)


def main():

    # import data as a list of dataframes
    raw_data = import_data(dat_path)

    # manipulate data
    cor_data = map(grass_correct_data, raw_data)
    new_data = map(find_ts_extrema, cor_data)
    ind_data = map(lambda x: get_extrema_points(x, tol=mytol), new_data)

    # now do optimization
    par_table = optimize_all_comb(ind_data)

    par_table.to_csv(out_path+"sigmoid_forcing.csv", index_label="k", columns=["value","error","site"])

if __name__ == '__main__':

    dat_path = "../data/"
    fig_path = "../figs/"
    out_path = "../outputs/"
    show_plot = False
    # Tolerance control on what points to extract around the extrema
    mytol = 1e-1

    main()


