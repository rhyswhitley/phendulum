#!/usr/bin/env python

import datetime, time
# load own modules
import spring_dynamics as _sd
import data_handling as _dh
import model_optim_extras as _mo
import model_plotting as _mp

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'


def main():
    dh = _dh.data_handling()
    mo = _mo.model_optim_extras()
    mp = _mp.model_plotting(fig_path)

    # import data as a list of Pandas dataframes
    raw_data = dh.import_data(dat_path)

    # map data transformations to each dataset imported
    cor_data = [ dh.grass_correct_data(rd) for rd in raw_data ]

    #p0 = [0.3, -10, 1, 1, 2]
    #bounds = [(0,10),(0,1000),(0,1),(-10,10),(-10,10)]

    p0 = [10, 0.1, 2, 10]
    bounds = [(0,1000),(0,1),(-10,10),(-10,10)]

    par_table = mo.optimize_all_sampling( mo.minimize_func, \
                    spring_motion, cor_data, p0, bounds, \
                    ylabel="NDVI_norm", xlabel="SWC10" )


    par_table.to_csv(out_path+"spring_parameters.csv", index_label="k",
                        columns=["Site","Sampling","Value","Error","EigVal"])

    par_casted = mo.recast_par_table(par_table)

    mp.plot_allSite_pendulum( cor_data, par_casted, mo.sig_mod1 )


def spring_motion(par, data):
    """ Function to return the prediction from the pendulum """
    mo = _mo.model_optim_extras()
    spring = _sd.spring(par, data, mo.sig_mod1, x_init=0.1, v_init=0.002 )
    motion = spring.calc_dynamics()['x']
    return motion

if __name__ == '__main__':

    dat_path = "../data/"
    fig_path = "../figs/"
    out_path = "../outputs/"
    tab_path = out_path + "parameter_table.csv"

    main()
else:
    print("Program FAIL: check config file")






#    def get_initials(vector):
#        x0 = vector[0]
#        v0 = vector[1] - x0
#        return (x0,v0)
#
#    xv_initials = [ get_initials(cd["NDVI_grass"]) for cd in cor_data]
#
#    init_spring = [ lambda par, data : \
#        _sd.spring(par, data, mo.sig_mod2, x_init=x0, v_init=v0) \
#        .calc_dynamics()['x']
#                   for (x0,v0) in xv_initials ]
#
#
