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
    ec_watr = ec_data.loc[:,("Precip_Con","Fn_Con","Fe_Con")]

    # Resample the hourly data to daily (although we could just do 16-day)
    phen_sampled = ec_phen.resample('D', how='mean',)
    phen_sampled.columns = ["SWC10","NDVI250X","VPD"]
    ndvi_pred = np.isfinite(phen_sampled['NDVI250X'])
    phen_filt = phen_sampled[ndvi_pred]

    temp_sampled = ec_temp.resample('D', how=('mean','min','max'),)
    temp_sampled.columns = ["Tmean","Tmin","Tmax"]
    temp_filt = temp_sampled[ndvi_pred]

    rain_sampled = ec_watr.resample('D', how='sum',)
    rain_sampled.columns = ["Rain","Rnet","AET"]
    rain_filt = rain_sampled[ndvi_pred]

    # net radiation cannot be less than 0
    rain_filt.Rnet[rain_filt.Rnet<0] = 1

    # put all preliminary information together here
    all_filt = pd.concat([phen_filt, temp_filt, rain_filt], axis=1)

    all_filt["Photoperiod"] = map( \
        lambda x: photoperiod(site_coord.Latitude,x), \
        phen_filt.index.dayofyear )

    all_filt["EET"] = [ equal_evap(T,A) for T,A in zip(all_filt.Tmean,all_filt.Rnet)]

    all_filt["Alpha"] = [ aet/eet for aet,eet in zip(all_filt.AET,all_filt.EET)]
    # correct for exceedingly high values
    all_filt.Alpha[all_filt.Alpha>1.26] = 1.26

    moo = np.array(all_filt["NDVI250X"])
    moo2 = f_frac(moo)

    plt.plot(all_filt["NDVI250X"], color='black', lw=2)
    plt.plot(moo2, color='red', lw=2)
    plt.show()

    #print all_filt.head()
    return None
    # Write to CSV into the Data folder
    #all_filt.to_csv(opath, sep=",")

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

def f_frac(ndvi):
    """
    4 stage process based on Donohue et al. (2009) to separate out tree and grass cover,
    using moving averages
    """
    # first calculate the 7-month moving minimum window across the time-series
    f1 = moving_something(min, ndvi)

    return f1

def moving_something(_fun, tseries, period=7, is_days=True):
    """
    Applies a function to a moving window of the time-series:
    ft_ = function([ f(t-N), f(t). f(t+N)])
    """
    # if the time-series is at a day-time step, update the window to a step-size of 16 days
    if is_days:
        p0 = period*16
    else:
        p0 = period

    # find upper and lower bounds of the moving window
    half = p0//2
    tlen = len(tseries)
    twin = [0]*tlen

    for t in tseries:
        # find the something for the window that satisfy the edge conditions
        if t < half:
            twin.append( _fun(tseries[t:t+half]) )
        elif t > tlen-half:
            twin.append( _fun(tseries[t-half:t]) )
        else:
            twin.append( _fun(tseries[t-half:t+half]) )

    return twin

def equal_evap(tair, Srad):
    s = sat_vap_slope(tair)
    g = psychometric(tair)
    return (s*Srad)/(s+g)

def sat_vap_slope(tair):
    """Slope of the relationships between vapour pressure and air temperature"""
    s = 6.1078*17.269*237.3*np.exp(17.269*tair/(237.3+tair))
    return 0.1*(s/(237.3+tair)**2)

def psychometric(tair):
    """Psychometric constant"""
    return 0.1*(0.646*np.exp(9.7E-4*tair))

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
#    natt_names = ["AdelaideRiver","AliceSprings","DalyUncleared","DryRiver", \
#                  "HowardSprings","SturtPlains"]
    natt_names = ["SturtPlains"]

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

