#!/usr/bin/env python

from os import listdir
import re, pickle, datetime, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
# own modules
import data_handling as _dh
import model_optim_extras as _mo
import spring_dynamics as _sd

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

def main():
    dh = _dh.data_handling()

    # import data as a list of Pandas dataframes
    raw_data = dh.import_data(dat_path)

    # map data transformations to each dataset imported
    cor_data = [ dh.normalise_xydata(rd) for rd in raw_data ]
    new_data = [ dh.find_ts_extrema(cd, var="NDVI_norm") for cd in cor_data ]
    ind_data = [ dh.get_extrema_points(nd, tol=0.1) for nd in new_data ]

    allofit = pd.concat(ind_data)
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.plot(allofit["SWC10"], allofit["NDVI_norm"], 'o')
    ax2.plot(allofit["RWC"], allofit["NDVI_norm"], 'o')
    plt.show()



    # import parameter posteriors
    post_table = pd.read_csv(out_path+ftable)

    return None

def plot_forcing(data):


    return None

if __name__=="__main__":

    out_path = "../outputs/pars/"
    dat_path = "../data/"
    ftable = "spring_posteriors.csv"
    main()
