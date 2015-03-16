#!/usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir
import re

pd.options.mode.chained_assignment = None

def main(fpath):

    # Display column names and their index for reference
    ec_data = pd.read_csv(fpath, parse_dates=True , index_col=['DT'])
    # import site coordinate information
    geocord = pd.read_csv(site_path)
    site_coord = geocord.loc[geocord["Label"]==site]

    if show_names == True:
        for n,i in zip(ec_data.columns.values,range(ec_data.shape[1])):
            print( str(i) + " : " + n )

    day_data = daily_agg(ec_data)
    all_data = calc_extras(day_data, site_coord)

    # Write to CSV into the Data folder
    all_data.to_csv(opath, sep=",", float_format='%11.6f')

def daily_agg(data):
    # Let just grab what we need -- SWC and NDVI
    ec_mean = data.loc[:,("Sws_Con","250m_16_days_NDVI_new_smooth","VPD_Con")]
    ec_range = data.loc[:,("Ta_Con")]
    ec_sums = data.loc[:,("Precip_Con","Fn_Con","Fe_Con")]

    # Resample the hourly data to daily (although we could just do 16-day)
    day_mean0 = ec_mean.resample('D', how='mean',)
    day_mean0.columns = ["SWC10","NDVI250X","VPD"]
    ndvi_pred = np.isfinite(day_mean0['NDVI250X'])
    day_mean1 = day_mean0[ndvi_pred]

    day_range0 = ec_range.resample('D', how=('mean','min','max'),)
    day_range0.columns = ["Tmean","Tmin","Tmax"]
    day_range1 = day_range0[ndvi_pred]

    day_sum0 = ec_sums.resample('D', how='sum',)
    day_sum0.columns = ["Rain","Rnet","AET"]
    day_sum1 = day_sum0[ndvi_pred]

    # net radiation cannot be less than 0
    day_sum1.loc[day_sum1['Rnet']<0, 'Rnet'] = 1

    # put all preliminary information together here
    all_filt = pd.concat([day_mean1, day_range1, day_sum1], axis=1)
    return all_filt

def calc_extras(data, geo):

    # determine the daily photoperiod
    data["Photoperiod"] = [ photoperiod(geo.Latitude,x) \
        for x in data.index.dayofyear ]

    # determine equilibrium evaporation
    data["EET"] = [ equal_evap(T,A) for T,A in zip(data.Tmean,data.Rnet)]

    # determine the Cramer-Prentice parameter 'Alpha'
    data["Alpha"] = [ aet/eet for aet,eet in zip(data.AET,data.EET)]
    # correct for exceedingly high values
    data.loc[data['Alpha']>1.26, 'Alpha'] = 1.26

    # determine the partitioning of NDVI between tree and grass layers
    tg_ratios = treegrass_frac(data["NDVI250X"], 30)
    data["NDVI_tree"] = tg_ratios["tree"]
    data["NDVI_grass"] = tg_ratios["grass"]

    return data

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

def treegrass_frac(ndvi, day_rs):
    """
    Process based on Donohue et al. (2009) to separate out tree and grass cover,
    using moving windows (adapted here for daily time-step)
    """
    # first calculate the 7-month moving minimum window across the time-series
    fp1 = moving_something(np.min, ndvi, period=7, day_rs=day_rs)
    fp2 = moving_something(lambda x: sum(x)/(9*day_rs), fp1, period=9, day_rs=day_rs)
    fr1 = ndvi - fp2

    ftree = [ p2-np.abs(r1) if r1<0 else p2 for p2,r1 in zip(fp2,fr1) ]
    fgrass = ndvi - ftree

    return ({'tree':ftree, 'grass':fgrass})

def moving_something(_fun, tseries, period, day_rs=16, is_days=True):
    """
    Applies a function to a moving window of the time-series:
    ft_ = function([ f(t-N), f(t). f(t+N)])
    """
    # if the time-series is at a day-time step, update the window to a step-size of 16 days
    if is_days:
        p0 = period*day_rs
    else:
        p0 = period

    # find upper and lower bounds of the moving window
    half = p0//2
    tlen = len(tseries)

    twin = [0]*tlen
    for i in range(tlen):
        # find the something for the window that satisfy the edge conditions
        if i < half:
            # fold back onto the end of the time-series
            twin[i] = _fun( np.hstack([ tseries[tlen-(half-i):tlen],\
                                        tseries[0:i+half]]) )
        elif i > tlen-half:
            # fold back into the beginning of the time-series
            twin[i] = _fun( np.hstack([ tseries[i-half:tlen],\
                                        tseries[0:half-(tlen-i)]]) )
        else:
            twin[i] = _fun(tseries[i-half:i+half])

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
    natt_names = ["AdelaideRiver","AliceSprings","DalyUncleared","DryRiver", \
                  "HowardSprings","SturtPlains"]
#    natt_names = ["SturtPlains"]

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

