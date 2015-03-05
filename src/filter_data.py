#!/usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir
import re

def main(fpath):

    # Display column names and their index for reference
    ec_data = pd.read_csv(fpath, parse_dates=True , index_col=['DT'])
    # import site coordinate information
    geocord = pd.read_csv(site_path)
    site_coord = geocord.loc[geocord["Label"]==site]

    if show_names == True:
        for n,i in zip(ec_data.columns.values,range(ec_data.shape[1])):
            print( str(i) + " : " + n )

    # Let just grab what we need -- SWC and NDVI
    ec_phen = ec_data.loc[:,("Sws_Con","250m_16_days_NDVI_new_smooth","VPD_Con")]
    ec_temp = ec_data.loc[:,("Ta_Con")]
    ec_rain = ec_data.loc[:,"Precip_Con"]

    # Resample the hourly data to daily (although we could just do 16-day)
    phen_sampled = ec_phen.resample('D', how='mean',)
    phen_sampled.columns = ["SWC10","NDVI250X","VPD"]
    ndvi_pred = np.isfinite(phen_sampled['NDVI250X'])
    phen_filt = phen_sampled[ndvi_pred]

    temp_sampled = ec_temp.resample('D', how=('mean','min','max'),)
    temp_sampled.columns = ["Tmean","Tmin","Tmax"]
    temp_filt = temp_sampled[ndvi_pred]

    rain_sampled = ec_rain.resample('D', how='sum',)
    rain_sampled.columns = ["Rainfall"]
    rain_filt = rain_sampled[ndvi_pred]

    # put it all together here
    all_filt = pd.concat([phen_filt, temp_filt, rain_filt], axis=1)

    all_filt["Photoperiod"] = map( \
        lambda x: photoperiod(site_coord.Latitude,x), \
        phen_filt.index.dayofyear )

    # Write to CSV into the Data folder
    all_filt.to_csv(opath, sep=",")

    # Add date/time index
    all_filt.reset_index(level=0, inplace=True)

    # Plots environmental data
    if draw_plot==True:

        fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True, figsize=(9,6))
        fig.subplots_adjust(hspace=0.1)
        fig.subplots_adjust(wspace=0.1)
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = "sans-serif"
        plt.rcParams['font.sans-serif'] = "Helvetica"
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['font.size'] = 12
        plt.rcParams['legend.fontsize'] = 9
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12

        almost_black = '#262626'
        # change the tick colors also to the almost black
        plt.rcParams['ytick.color'] = almost_black
        plt.rcParams['xtick.color'] = almost_black

        # change the text colors also to the almost black
        plt.rcParams['text.color'] = almost_black

        # Change the default axis colors from black to a slightly lighter black,
        # and a little thinner (0.5 instead of 1)
        plt.rcParams['axes.edgecolor'] = almost_black
        plt.rcParams['axes.labelcolor'] = almost_black

        ax1.plot(all_filt['DT'], all_filt["NDVI250X"], linestyle='-', color='red')
        ax2.plot(all_filt['DT'], all_filt["SWC10"], linestyle='-', color='blue')
        ax1.set_ylabel(r"NDVI (-)")
        ax2.set_ylabel(r"$\theta_{10cm}$")
        ax2.set_xlabel(r"Years")

        # force zero start in range
        ax1.set_ylim(ymin=0)
        ax2.set_ylim(ymin=0)

        ax1.yaxis.major.locator.set_params(nbins=5)
        ax2.yaxis.major.locator.set_params(nbins=5)

        simpleaxis(ax1)
        simpleaxis(ax2)

        #fig.autofmt_xdate()
        fig.savefig("{0}{1}_filt.pdf".format(fig_path,site),
                    bbox_inches='tight', pad_inches=0.1)

def simpleaxis(ax):
    """ Remove the top line and right line on the plot face """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def get_all_site_names(ec_files):
    # list comprehensive way of getting all names in bulk
    file_name = [ re.compile('^\w+').findall(f) for f in ec_files ]
    split_name = [ re.split('_',f) for f in sum(file_name,[]) ]
    name_ver = [[ name for i,name in enumerate(sublist) if (i==3) ] \
                  for sublist in split_name ]
    return sum(name_ver,[])

def get_site_name(ec_files):
    file_name = re.compile('^\w+').findall(ec_files)
    split_name = re.split('_',file_name[0])
    return split_name[-2]

def photoperiod(lat, doy):
    """Ecological Modeling, volume 80 (1995) pp. 87-95"""
    P = math.asin(.39795*math.cos(.2163108 + \
                    2*math.atan(.9671396*math.tan(.00860*(doy-186)))))
    numerator = math.sin(0.8333*math.pi/180.) + \
                    math.sin(lat*math.pi/180.)*math.sin(P)
    denominator = math.cos(lat*math.pi/180.)*math.cos(P)
    photohr = 24. - (24./math.pi)*math.acos(numerator/denominator)
    return photohr


if __name__ == '__main__':

    # Import EC tower Dingo dataset for Sturt Plains [move this to a config file]
    out_path = "../data/"
    fig_path = "../figs/"
    search_path = out_path+"Dingo_v12/"
    site_path = search_path+"site coordinates ozflux.csv"
    show_names = False
    draw_plot = True
    out_name = "filtered"
    version = "_v12"

    # collect all files for processed eddy covariance datasets in the data folder
    natt_names = ["AdelaideRiver","AliceSprings","DalyUncleared","DryRiver", \
                  "HowardSprings","SturtPlains"]

    get_files = [ f for f in listdir(search_path) if f.endswith('.csv') ]

    grab_these = pd.Series(get_files).str.contains("|".join(natt_names))

    natt_files = pd.Series(get_files)[np.array(grab_these)]

    # for each site extract the required data
    for i_file in natt_files:
        site = get_site_name(i_file)
        print("Filtering data for site => "+ site)
        opath = out_path+out_name+"_"+site+version+".csv"
        fpath = search_path+i_file
        main(fpath)

