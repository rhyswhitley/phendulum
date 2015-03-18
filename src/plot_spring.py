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

    # import parameter posteriors
    post_table = pd.read_csv(out_path+ftable)

    # create plots
    no_figs = len(cor_data)
    fig = plt.figure(figsize=(12,9))
    gs = gridspec.GridSpec(3,2)

    sax = fig.add_subplot(111) # <<< big subplot
    # Turn off axis lines and ticks of the big subplot
    sax.spines['top'].set_color('none')
    sax.spines['bottom'].set_color('none')
    sax.spines['left'].set_color('none')
    sax.spines['right'].set_color('none')
    sax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    for i,data in enumerate(cor_data):
        pObj = fig.add_subplot(gs[i])
        stab = post_table[post_table.Site==data["Site"][0]]
        if len(stab)<5:
            # get initial conditions
            x_init = data.NDVI_norm[0]
            v_init = data.NDVI_norm[1] - data.NDVI_norm[0]
            new_mean = np.concatenate([stab.Mean,[x_init,v_init]], axis=0)
            new_ci05 = np.concatenate([stab.CI05,[x_init,v_init]], axis=0)
            new_ci95 = np.concatenate([stab.CI95,[x_init,v_init]], axis=0)
            _plotSpring(pObj, data, new_mean, new_ci05, new_ci95, stab.Site.iloc[0])
        else:
            _plotSpring(pObj, data, stab.Mean, stab.CI05, stab.CI95, stab.Site.iloc[0])
    sax.set_ylabel("Approximated Grass NDVI Response", size=14)
    sax.yaxis.labelpad = 20
    plt.tight_layout(h_pad=2)
    plt.savefig("../figs/phendulum_bayes_exp.pdf")


def _plotSpring(ax, data, k_means, k_05, k_95, site_name, xlabel="SWC10"):

    years = YearLocator()   # every year
    months = MonthLocator()  # every month
    yearsFmt = DateFormatter('%Y')

    mo = _mo.model_optim_extras()
    eforce = mo.sig_mod1

    prediction = _sd.spring(k_means, data[xlabel], eforce) \
                    .calc_dynamics()['x']
    pred_ci_U = _sd.spring( k_95, data[xlabel], eforce) \
                    .calc_dynamics()['x']
    pred_ci_L = _sd.spring( k_05, data[xlabel], eforce) \
                    .calc_dynamics()['x']

    ax.fill_between( data.index, pred_ci_U,  pred_ci_L, color='red', alpha=0.2)
    ax.plot_date( data.index[::16], data["NDVI_norm"][::16], 'o', linestyle='-', lw=1.2, color="black", alpha=0.6)
    ax.plot_date( data.index, prediction, '-', lw=2.5, color="#DC143C")

    ax.set_title(site_name, size=16)

    # format the ticks
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(months)
    ax.autoscale_view()

    plt.setp(plt.xticks()[1], rotation=30, ha='right')

if __name__=="__main__":

    out_path = "../outputs/pars/"
    dat_path = "../data/"
    ftable = "spring_posteriors.csv"
    main()


