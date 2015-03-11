#!/usr/bin/env python

from os import listdir
import re, pickle, datetime, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
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
    cor_data = [ dh.grass_correct_data(rd) for rd in raw_data ]

    # import parameter posteriors
    post_table = pd.read_csv(out_path+ftable)

    # create plots
    for data in cor_data:
        #print post_table["Site"==data["Site"][0]]
        stab = post_table[post_table.Site==data["Site"][0]]
        print pd.concat(stab.Mean,[0,0])
        #_plotSpring(data, stab.Mean, stab.CI05, stab.CI95)


def _plotSpring(data, k_means, k_05, k_95, xlabel="SWC10"):

    mo = _mo.model_optim_extras()
    eforce = mo.sig_mod1

    prediction = _sd.spring(k_means, data[xlabel], eforce) \
                    .calc_dynamics()['x']
    pred_ci_U = _sd.spring( k_95, data[xlabel], eforce) \
                    .calc_dynamics()['x']
    pred_ci_L = _sd.spring( k_05, data[xlabel], eforce) \
                    .calc_dynamics()['x']

    plt.fill_between( range(len(prediction)), pred_ci_U,  pred_ci_L, color='red', alpha=0.2)
    plt.plot( data["NDVI_norm"], '.', lw=2, color="black")
    plt.plot( prediction, lw=2, color="#DC143C")
    plt.axis([0,len(prediction),0,1])
    plt.show()

if __name__=="__main__":

    out_path = "../outputs/"
    dat_path = "../data/"
    ftable = "spring_posteriors.csv"
    main()


