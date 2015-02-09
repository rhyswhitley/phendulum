#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import datetime, time
# import our own modules
import data_handling as _dh
import springDynamics as _sd

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'


def opt_environmental_forcing(raw_data):
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

    # now do an optimization on all the imported datasets
    par_table = dh.optimize_all_sites(ind_data)

    # create a comma delimited table of optimised environmental forcing
    par_table.to_csv(out_path+"sigmoid_forcing.csv", index_label="k", \
                     columns=["value","error","site"])

def opt_phendulum(raw_data):


# Main
def main():
    # import data as a list of Pandas dataframes
    raw_data = _dh.data_handling.import_data(dat_path)

    # creates the parameter table that describes the environmental forcing used on the spring
    opt_environmental_forcing(raw_data)

if __name__ == '__main__':
    dat_path = "../data/"
    fig_path = "../figs/"
    out_path = "../outputs/"
    show_plot = False
    # Tolerance control on what points to extract around the extrema
    mytol = 1e-1
    main()
else:
    print "FAIL"







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

