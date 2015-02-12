#!/usr/bin/env python

import datetime, time
import pandas as pd
# load own modules
import spring_dynamics as _sd
import data_handling as _dh
import model_optim_extras as _mo
#import model_plotting as _mp

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'


def main():
    dh = _dh.data_handling()
    mo = _mo.model_optim_extras()

    # import data as a list of Pandas dataframes
    raw_data = dh.import_data(dat_path)

    # map data transformations to each dataset imported
    cor_data = [ dh.grass_correct_data(rd) for rd in raw_data ]

    par_table = [ mo.optimize_func( \
                    spring_motion, _data, p0=[-10,-3,1,1], \
                    ylabel="NDVI_grass", xlabel="SWC10" )
                 for _data in cor_data ]

#    par_table = mo.optimize_all_sampling(
#                    spring_motion, cor_data, p0=[-10,-3,1,1], \
#                    ylabel="NDVI_grass", xlabel="SWC10" )

    print pd.concat(par_table)

def spring_motion(par, data):
    mo = _mo.model_optim_extras()
    spring = _sd.spring(par, data, mo.sig_mod1)
    motion = spring.calc_dynamics()['x']
    return motion

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
    tb_pivot = map( lambda x: x.reset_index().pivot(index="Label", columns="k",
                    values="Value"), tab_sub )
    # return to user
    return tb_pivot


if __name__ == '__main__':

    dat_path = "../data/"
    fig_path = "../figs/"
    out_path = "../outputs/"
    tab_path = out_path + "parameter_table.csv"

    main()
else:
    print "Program FAIL: check config file"







