#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime, time
# load own modules
import data_handling as _dh
import model_optim_extras as _mo
import springDynamics as _sd

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'


def opt_environmental_forcing(raw_data, find_params=False):
    """
    Determines the environmental forcing that drives the momentum of the system
    (or pendulum). Return a table of fitted parameter values that describe the
    shape of the forcing based on the local environmental data.
    """
    # create a data handler object
    dh = _dh.data_handling()

    # map data transformations to each dataset imported
    cor_data = map(dh.grass_correct_data, raw_data)
    new_data = map(dh.find_ts_extrema, cor_data)
    ind_data = map(lambda x: dh.get_extrema_points(x, tol=mytol), new_data)

    if find_params:
        # now do an optimization on all the imported datasets
        par_table = dh.optimize_on_sampling(ind_data)

        # create a comma delimited table of optimised environmental forcing
        par_table.to_csv(out_path+"sigmoid_forcing.csv", index_label="k", \
                        columns=["Value","Error","Site","Sampling"])
        return [new_data, par_table]
    else:
        # import parameters table from file
        file_table = pd.read_csv(tab_path)
        return [new_data, file_table]

def recast_par_table(tab_im):
    """
    Function accepts a parameter table as derived by the environmental
    forcing procedure and re-casts as a list of parameter values for each
    sampling type
    """
    # get the prefixes for each site name (stored in the ensemble site name)
    ens_lab = tab_im.query('Sampling=="ensemble"').Site[0]
    prefix = ens_lab.split("_")
    # add new column to pivot on
    tab_im["Label"] = tab_im["Site"]+"_"+tab_im["Sampling"]
    # now need to group the parameter table into site subsets
    tab_sub = map( lambda x: tab_im[tab_im["Site"].str.contains(x)], prefix )
    # recast the above list of subsets into a list of callable parameters
    tb_pivot = map( lambda x: x.pivot(index="Label", columns="k", values="Value"), tab_sub )
    # return to user
    return tb_pivot

# Main routine
def main():

    mo = _mo.model_optim_extras()

    # import data as a list of Pandas dataframes
    raw_data = _dh.data_handling().import_data(dat_path)

    # creates the parameter table that describes the environmental forcing used on the spring
    data_list, raw_table = opt_environmental_forcing(raw_data, find_params=False)

    # turns the flat parameter table into N-D list of site parameter values at different samplings
    par_cast = recast_par_table(raw_table)
    # make sure that the number of grouped parameters in the dataset equals the number of imported datasets
    assert len(data_list)==len(par_cast), "Number of datasets exceeds the number of available parameter sets"

    # transform list of pivot tables into a list of parameter list for each site and sampling
    par_list = map( lambda x: map(list,np.array(x)), par_cast)
    # now find expressions for the environmental forcing at each site for each sampling type
    force_list = [ [mo.environ_force(k,mo.sig_mod) for k in k_val] for k_val in par_list ]
    # create a list of springs based on the above list of forces
    springs_list = [ [create_spring(f) for f in fs] for fs in force_list ]
    # now fit these springs to the NDVI time-series for each site
    spring_coeff = [ [mo.optimise_func( sub_spring, data, p0=[1,1]) for sub_spring in springs] for (springs,data) in zip(springs_list,data_list) ]
    print spring_coeff

def create_spring(force_on):
    fx_model = lambda par, X: _sd.spring( par, X, force_on ).calc_dynamics()['x']
    return fx_model


if __name__ == '__main__':
    dat_path = "../data/"
    fig_path = "../figs/"
    out_path = "../outputs/"
    tab_path = out_path + "sigmoid_forcing.csv"
    show_plot = False
    # Tolerance control on what points to extract around the extrema
    mytol = 1e-1
    main()
else:
    print "Program FAIL: check config file"







def plot_force_optfits(dataset, p_table, fe_mod, \
                        xlabel="SWC_smooth", ylabel="NDVI_grass", \
                        file_name = "_phen_fe_fit.pdf"):
    # Create vectors for model fits
    xs = np.arange(0,0.3,1e-3)
#        ndvi_lin = lin_mod( p_table[p_table["lin"]]['value'], xs )
#        ndvi_exp = exp_mod( p_table[p_table["exp"]]['value'], xs )
    ndvi_sig = fe_mod( p_table[p_table["sig"]]['value'], xs )
    # needs to be a list comprehension
    ndvi_sig

    # Plot the results
    plt.plot( dataset[xlabel], dataset[ylabel], 'o', color='black' )
    plt.plot( xs, ndvi_lin, linestyle='-', color='red', lw=2, label=r"$k_{0}\theta_{s}$" )
    plt.plot( xs, ndvi_exp, linestyle='-', color='blue', lw=2, label=r"$k_{1}\exp(k_{2}\theta_{s})$" )
    plt.plot( xs, ndvi_sig, linestyle='-', color='purple', lw=2, label=r"$(1+\exp(k_{3}\theta_{s}-k_{4}))^{-1}$" )
    plt.xlabel(r'$\theta_{s 10cm}$', fontsize=18)
    plt.ylabel('NDVI')
    plt.legend(loc=2)
    #plt.axis([0,0.25,0,0.5])
    #self.is_plotted(fig, file_name)
    plt.show()

