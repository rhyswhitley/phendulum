#!/usr/bin/env python

import numpy as np
import pandas as pd
import re, math
import datetime, time
from os import listdir
from scipy.ndimage import gaussian_filter

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015, 1, 14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

class data_handling(object):
    def __init__(self):
        None
    def import_data(self, dat_path, includes=r'^filtered.*'):
        """
        Function reads in datasets stored in the /data folder in the upper
        directory, attaches site names based on a regex search on the file names
        and concatenates the list of datasets into one Pandas dataframe for
        analysis
        """
        # get file names
        file_names = [ f for f in listdir(dat_path) if re.match(includes,f) ]
        assert file_names!=[], "There are no filtered datasets in the "+dat_path+" location: run the filter script first"
        # read files in
        datasets = [ pd.read_csv(dat_path+x, parse_dates=True, index_col=['DT']) for x in file_names ]
        # find site names
        site_names = self.get_all_site_names( file_names )
        # attach site names relating to each imported dataset
        for i, df in enumerate(datasets):
            df['Site'] = site_names[i]
        # concatenate and return final dataframe
        return datasets

    def normalise_xydata(self, dataset, xvar="SWC10", yvar="NDVI250X"):
        """
        Function filters data to enable a better approximation of the
        environmental forcing that drives the momentum of the pendulum. In this
        case we have smoothed the soil moisture response to reduce the number of
        extrema points and better align with the extrema in the NDVI time-series
        """
        x_raw = dataset[xvar]
        y_raw = dataset[yvar]
        # smooth SMC to obtain more identifiable extrema points (removes noise)
        x_approx = gaussian_filter(x_raw, 10)
        # normalise the NDVI grass estimate to bring into line with Jolly
        y_norm = (y_raw - min(y_raw))/(max(y_raw) - min(y_raw))
        # normalise the SWC to account for differences in soil chars across sites
        x_norm = (x_raw - min(x_raw))/(max(x_raw) - min(x_raw))

        # assign new values as extra columns to the dataframe
        dataset['SWC_smooth'] = x_approx
        dataset['RWC'] = x_norm
        dataset['NDVI_norm'] = y_norm
        # return to user
        return dataset

    def find_ts_extrema(self, dataset, var="NDVI_grass"):
        """
        Finds the maxima and minima or turning points of a time-series
        """
        # Find the inflexion points in the NDVI time-series
        dataset['dy1/dt1'] = dataset[var].shift(1).fillna(0)-dataset[var]
        dataset['dy2/dt2'] = dataset['dy1/dt1'].shift(1).fillna(0)-dataset['dy1/dt1']
        return dataset

    def get_extrema_points(self, dataset, tol=1):
        """
        Extract data points along the dy/dx time-series that intersects with the
        dy2/dx2 time-series. This should approximate points that occur around the
        extrema of the dependent variable time-series
        """
        # Put these in a dataframe to filter sp_data_new
        upper = max(dataset["dy2/dt2"].ix[2:])*tol
        lower = min(dataset["dy2/dt2"].ix[2:])*tol
        dataset_ex = dataset.loc[(dataset['dy1/dt1']<upper) & (dataset['dy1/dt1']>lower)]
        if len(dataset_ex.index) <= 1:
            print("ERROR: Tolerance is too high :: not enough data to optimize on")
            return None
        else:
            return dataset_ex #.ix[:,'SWC_smooth':'NDVI250X']

    def get_all_site_names(self, ec_files, regex='^\w+', pos=2):
        """
        Extract the site names from an array of file names assuming they follow the
        Dingo file name format, i.e. Advanced_processed_SITENAME_vXX.csv
        """
        file_name = (re.compile(regex).findall(f) for f in ec_files)
        split_name = (re.split('_', f) for f in sum(file_name, []))
        name_ver = [sublist[-pos] for sublist in split_name]
        return name_ver

    def df_pop_site(self, dseries):
        """
        Return unique key values from dictionary
        """
        return [ s for s in set(dseries) ]

    def create_label(self, dseries):
        """
        If there is more than one unique label value, then take the first two
        characters and join them into a new string with an underscore linker.
        Otherwise just return the singular name as the label.
        """
        sites = self.df_pop_site(dseries)
        if len(sites) > 1:
            prefix = [char[:2] for char in sites]
            return '_'.join(prefix)
        else:
            return sites.pop()

    def get_all_site_names(self, ec_files):
        # list comprehensive way of getting all names in bulk
        file_name = [re.compile('^\w+').findall(f) for f in ec_files]
        split_name = [re.split('_',f) for f in sum(file_name,[])]
        name_ver = [[name for i,name in enumerate(sublist) if (i==3)] \
                    for sublist in split_name]
        return sum(name_ver,[])

    def get_site_name(self, ec_files):
        file_name = re.compile('^\w+').findall(ec_files)
        split_name = re.split('_', file_name[0])
        return split_name[-2]

    def photoperiod(self, lat, doy):
        """Ecological Modeling, volume 80 (1995) pp. 87-95"""
        P = math.asin(.39795*math.cos(.2163108 + \
                        2*math.atan(.9671396*math.tan(.00860*(doy-186)))))
        numerator = math.sin(0.8333*math.pi/180.) + \
                        math.sin(lat*math.pi/180.)*math.sin(P)
        denominator = math.cos(lat*math.pi/180.)*math.cos(P)
        photohr = 24. - (24./math.pi)*math.acos(numerator/denominator)
        return photohr

    def treegrass_frac(self, ndvi, day_rs):
        """
        Process based on Donohue et al. (2009) to separate out tree and grass cover,
        using moving windows (adapted here for daily time-step)
        """
        # first calculate the 7-month moving minimum window across the time-series
        fp1 = self.moving_something(np.min, ndvi, period=7, day_rs=day_rs)
        fp2 = self.moving_something(lambda x: sum(x)/(9*day_rs), fp1, period=9, day_rs=day_rs)
        fr1 = ndvi - fp2

        ftree = [p2 - np.abs(r1) if r1 < 0 else p2 for (p2, r1) in zip(fp2, fr1)]
        fgrass = ndvi - ftree

        return ({'tree':ftree, 'grass':fgrass})

    def moving_something(self, _fun, tseries, period, day_rs=16, is_days=True):
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
                twin[i] = _fun(np.hstack([tseries[tlen-(half-i):tlen],\
                                            tseries[0:i+half]]))
            elif i > tlen-half:
                # fold back into the beginning of the time-series
                twin[i] = _fun(np.hstack([tseries[i-half:tlen],\
                                            tseries[0:half-(tlen-i)]]))
            else:
                twin[i] = _fun(tseries[i-half:i+half])

        return twin

    def equal_evap(self, Tair, Srad):
        """Determination of the available energy or equilibrium evaporation"""
        # slope of the relationship between temperature and vapour pressure deficit
        slope = self.sat_vap_slope(Tair)
        # pyschometric constant
        psych = self.psychometric(Tair)
        # calculate available energy
        return (slope*Srad)/(slope + psych)

    def sat_vap_slope(self, tair):
        """Slope of the relationships between vapour pressure and air temperature"""
        s = 6.1078*17.269*237.3*np.exp(17.269*tair/(237.3+tair))
        return 0.1*(s/(237.3+tair)**2)

    def psychometric(self, tair):
        """Psychometric constant"""
        return 0.1*(0.646*np.exp(9.7E-4*tair))

