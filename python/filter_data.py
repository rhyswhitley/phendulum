#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os.path import expanduser, isfile

# Import EC tower Dingo dataset for Sturt Plains
out_fold = "data/"
out_name = "filtered"
site = "SturtPlains"
version = "_v12"
opath = "{0}{1}_{2}.csv".format(out_fold,out_name,site)
fpath = "data/Advanced_processed_data_{0}{1}.csv".format(site,version)

if isfile(opath)==True:
    print("File already exists")
else:
    # Display column names and their index for reference
    ec_data = pd.read_csv( fpath, parse_dates=True , index_col=['DT'] )
    show_names = False
    if show_names == True:
        for n,i in zip(ec_data.columns.values,range(ec_data.shape[1])):
            print( str(i) + " : " + n )

    # Let just grab what we need -- SWC and NDVI
    ec_phen = ec_data.loc[:,("Sws_Con","250m_16_days_NDVI_new_smooth")]

    # Resample the hourly data to daily (although we could just do 16-day)
    ec_filt = ec_phen.resample('D', how='mean',)
    ec_filt.columns = ["SWC10","NDVI250X"]

    show_plot = True
    # Plots environmental data
    if show_plot==True:
        fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot( ec_filt["NDVI250X"], linestyle='-', color='red' )
        ax2.plot( ec_filt["SWC10"], linestyle='-', color='blue' )
        ax1.set_ylabel( r"NDVI", fontsize=14 )
        ax2.set_ylabel( r"$\theta_{s 10cm}$", fontsize=18 )
        #fig.autofmt_xdate()
        plt.show()

    # Write to CSV into the Data folder
    ec_filt.to_csv( opath, sep="," )
