#!/usr/bin/env python

import re, os
import datetime, time
import pandas as pd
from scipy import integrate

from data_handling import data_handling

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015, 1, 14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

pd.options.mode.chained_assignment = None

def daily_agg(data):
    """
    Cleans and aggregates eddy-covariance data to a daily frequency that
    matches the NDVI observations
    """

    # lambda function for integrated sum of daily energy
    dayint_MJ = lambda x: integrate.trapz(x, dx=1800)*1e-6

    # wanted variables from full set
    ec_vars = ["Sws_Con", "250m_16_days_NDVI_new_smooth", "VPD_Con", \
                "Ta_Con", "Precip_Con", "Fn_Con", "Fe_Con"]

    # only take the columns we need
    sample_set = data[ec_vars]
    # get rid of NaNs
    sample_set.fillna(method='bfill', inplace=True)
    # clean the dataset for extraneous radiation values
    sample_set.ix[sample_set['Fn_Con'] < 0.0, ['Fn_Con']] = 0.0
    # rename columns (easier to work with)
    sample_set.columns = ["SWC10", "NDVI250X", "VPD", "Tair", "Rain", "Rnet", "AET"]
    # add equilibrium evaporation
    sample_set["EET"] = [DH.equal_evap(Ta, Rn) if Rn > 0 else 1e-6 for (Ta, Rn) \
                    in zip(sample_set.Tair, sample_set.Rnet)]

    # create a dictionary to resample these variables at a daily frequency
    sampling = {lab: 'mean' if i < 3 else ('min', 'max', 'mean') \
                if i == 3 else 'sum' if i == 4 else dayint_MJ \
                for (i, lab) in enumerate(sample_set.columns)}

    # re-sample data to a daily time-series
    daily_tseries = sample_set.resample('D', how=sampling)
    # drop the multi-index columns (but keep the multiple Tair labels)
    daily_tseries.columns = ["{0}_{1}".format(p1, p2) \
                            if p1 == 'Tair' else p1 \
                            for (p1, p2) in daily_tseries.columns.values]
    return daily_tseries

def calc_extras(data, geo):
    """
    Extra calculations that are necessary to operate the phenology
    pendulum.

    Photoperiod: is not used yet, but may be useful for northern hemisphere
    studies.

    Alpha: acts as a proxy for soil moisture content in the upper soil surface.
    Required for applying the phendulum spatially.

    NDVI_(grass/tree): Determines the recurrent and persistent signals of the
    bulk LAI signal. Currently the phendulum only operates to predict the
    life-cycle of the recurrent signal.
    """

    # determine the daily photoperiod
    data["Photoperiod"] = [DH.photoperiod(geo.Latitude, x) \
                    for x in data.index.dayofyear]

    # determine the Cramer-Prentice parameter 'Alpha'
    data["Alpha"] = [max(0, min(1.26, ET/Eq)) \
                    for (ET, Eq) in zip(data.AET, data.EET)]

    # determine the partitioning of NDVI between tree and grass layers
    tg_ratios = DH.treegrass_frac(data["NDVI250X"], 30)
    data["NDVI_tree"] = tg_ratios["tree"]
    data["NDVI_grass"] = tg_ratios["grass"]

    return data

def main(fpath):
    """
    Script that reduces any half-hourly OzFlux datasets to a daily aggregated
    version that contains the information needed to calibrate phenology models
    that describe bioclimate thresholds for leaf emergence and abscission. I.e
    the datasets produced from this script should ideally be capable of
    calibrating both Jolly et al (2005) and the phendulum model.

    This of course assumes that all OzFlux datasets adhere to the same naming
    convention for each *measured* variable - something I doubt.

    MAIN processes are run below, and make sure the data_handling folder (or
    symlink to it) exists in the current directory.
    """

    # Display column names and their index for reference
    ec_data = pd.read_csv(fpath, parse_dates=True, index_col=['DT'])

    # import OzFlux site geospatial information (needed for photoperiod)
    geocord = pd.read_csv(SITE_PATH)

    # only get the site coordinates for the one being processed
    site_coord = geocord.loc[geocord["Label"] == SITE]

    # determine the daily time-series from the half-hourly data
    day_data = daily_agg(ec_data)

    # determine extra variables (tree/grass NDVI, Alpha, photoperiod)
    all_data = calc_extras(day_data, site_coord)

    # Write to CSV into the Data folder
    all_data.to_csv(CSVPATH, sep=",", float_format='%11.6f')

    return None

if __name__ == '__main__':

    # create a data handling object
    DH = data_handling()

    # Import EC tower Dingo dataset for Sturt Plains [move this to a config file]
    OUT_PATH = "../data/"
    SEARCH_PATH = OUT_PATH + "Dingo_v12/"
    SITE_PATH = SEARCH_PATH + "site coordinates ozflux.csv"
    OUT_NAME = "filtered"
    VERSION = "_v12"

    # collect all files for processed eddy covariance datasets in the data folder
    NATT_NAMES = ["AdelaideRiver", "AliceSprings", "DalyUncleared", "DryRiver", \
                  "HowardSprings", "SturtPlains"]

    # only get CSV files
    NATT_FILES = [os.path.join(dp, f) for (dp, _, fn) in os.walk(SEARCH_PATH) \
                  for f in fn if re.search("|".join(NATT_NAMES), f)]

    # for each site extract the required data
    for nfile in NATT_FILES:

        SITE = DH.get_site_name(os.path.basename(nfile))
        print("Filtering data for site => " + SITE)

        CSVPATH = OUT_PATH + OUT_NAME + "_" + SITE + VERSION + ".csv"

        main(nfile)

