#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os.path import isfile
import sys

# Import EC tower Dingo dataset for Sturt Plains [move this to a config file]
out_fold = "../data/"
fig_fold = "../figs/"
out_name = "filtered"
site = "SturtPlains"
version = "_v12"
opath = "{0}{1}_{2}{3}.csv".format(out_fold,out_name,site,version)
fpath = "{0}{1}{2}.csv".format(out_fold,site,version)

if isfile(opath)==True:
    print("File already exists")
else:
    # Display column names and their index for reference
    ec_data = pd.read_csv( fpath, parse_dates=True , index_col=['DT'] )
    show_names = False
    show_plot = True
    if show_names == True:
        for n,i in zip(ec_data.columns.values,range(ec_data.shape[1])):
            print( str(i) + " : " + n )

    # Let just grab what we need -- SWC and NDVI
    ec_phen = ec_data.loc[:,("Sws_Con","250m_16_days_NDVI_new_smooth")]

    # Resample the hourly data to daily (although we could just do 16-day)
    ec_sampled = ec_phen.resample('D', how='mean',)
    ec_sampled.columns = ["SWC10","NDVI250X"]
    ec_filt = ec_sampled[np.isfinite(ec_sampled['NDVI250X'])]

    draw_plot = True
    # Plots environmental data
    if draw_plot==True:
        fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot( ec_filt["NDVI250X"], linestyle='-', color='red' )
        ax2.plot( ec_filt["SWC10"], linestyle='-', color='blue' )
        ax1.set_ylabel( r"NDVI", fontsize=14 )
        ax2.set_ylabel( r"$\theta_{s 10cm}$", fontsize=18 )
        #fig.autofmt_xdate()
        fig.savefig("{0}{1}_filt.pdf".format(fig_fold,site))

    # Write to CSV into the Data folder
    ec_filt.to_csv( opath, sep="," )
