#!/usr/bin/env python

import datetime, time
import pandas as pd
# load own modules
import data_handling as _dh
import model_optim_extras as _mo
import spring_dynamics as _sd
import model_plotting as _mp

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'


def main():


    # import data as a list of Pandas dataframes
    raw_data = _dh.data_handling().import_data(dat_path)

    # creates the parameter table that describes the environmental forcing used
    # on the spring
    (data_list, extd_list,
     raw_table) = opt_environmental_forcing(e_force, raw_data, find_params=True)

    # turns the flat parameter table into N-D list of site parameter values at
    # different samplings
    par_cast = mo.recast_par_table(raw_table)

    # make sure that the number of grouped parameters in the dataset equals the
    # number of imported datasets
    msg = "Number of datasets exceeds the number of available parameter sets"
    assert len(data_list)==len(par_cast), msg

    mplots = _mp.model_plotting(fig_path)
    mplots.plot_data_manipulation(data_list)
    mplots.plot_allSite_forcing(e_force, extd_list, par_cast)

def opt_environmental_forcing(f_mod, raw_data, find_params=False):
    """
    Determines the environmental forcing that drives the momentum of the system
    (or pendulum). Return a table of fitted parameter values that describe the
    shape of the forcing based on the local environmental data.
    """
    # create objects
    dh = _dh.data_handling()
    mo = _mo.model_optim_extras()

    # map data transformations to each dataset imported
    cor_data = [ dh.grass_correct_data(rd) for rd in raw_data ]
    new_data = [ dh.find_ts_extrema(cd) for cd in cor_data ]
    ind_data = [ dh.get_extrema_points(nd, tol=mytol) for nd in new_data ]

    if find_params:
        # now do an optimization on all the imported datasets
        bounds = [(0,1),(0,1000),(0,1)]
        p0 = [0.3,10,0.1]
        par_table = mo.optimize_all_sampling(mo.optimize_func,
                                             f_mod, ind_data, p0, bounds,
                                             ylabel="NDVI_grass",
                                             xlabel="SWC_smooth")

        # create a comma delimited table of optimised environmental forcing
        par_table.to_csv(out_path+"sigmoid_forcing.csv", index_label="k",
                         columns=["Value","Error","Site","Sampling"])
        return [new_data, ind_data, par_table]
    else:
        # import parameters table from file
        try:
            file_table = pd.read_csv(tab_path, index_col=0)
            return [new_data, ind_data, file_table]
        except IOError:
            print "\nNo parameter file exists: will create one now\n"
            # call recursively; find_params=True is the edge condition
            return opt_environmental_forcing(raw_data, f_mod, find_params=True)

if __name__ == '__main__':

    dat_path = "../data/"
    fig_path = "../figs/"
    out_path = "../outputs/"
    tab_path = out_path + "sigmoid_forcing.csv"
    show_plot = False
    # set the type of external forcing model here
    mo = _mo.model_optim_extras()
    e_force = mo.sig_mod2

    # Tolerance control on what points to extract around the extrema
    mytol = 5e-2
    main()
else:
    print "Program FAIL: check config file"






