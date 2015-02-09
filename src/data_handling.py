#!/usr/bin/env python

import pandas as pd
import re
import datetime, time
from os import listdir
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter
# import our own modules
import model_optim_extras as _mo

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

class data_handling(object):

    #def __init__(self, dat_path=""):
    #    self.dpath = dat_path

    def import_data(self, dat_path, includes=r'^filtered.*'):
        """
        Function reads in datasets stored in the /data folder in the upper
        directory, attaches site names based on a regex search on the file names
        and concatenates the list of datasets into one Pandas dataframe for
        analysis
        """
        # get file names
        file_names = [ f for f in listdir(dat_path) if re.match(includes,f) ]
        # read files in
        datasets = [ pd.read_csv(dat_path+x) for x in file_names ]
        # find site names
        site_names = self.get_all_site_names( file_names )
        # attach site names relating to each imported dataset
        for i, df in enumerate(datasets):
            df['Site'] = site_names[i]
        # concatenate and return final dataframe
        return datasets

    def grass_correct_data(self, dataset, xvar="SWC10", yvar="NDVI250X"):
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
        # adjust NDVI by subtracting minimum time-series value to approximate grasses
        y_approx = y_raw - min(y_raw)
        # assign new values as extra columns to the dataframe
        dataset['SWC_smooth'] = x_approx
        dataset['NDVI_grass'] = y_approx
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

    def get_extrema_points(self, dataset, tol=1e-2):
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
            print "ERROR: Tolerance is too high :: not enough data to optimize on"
            return None
        else:
            return dataset_ex

    def get_all_site_names(self, ec_files):
        """
        Extract the site names from an array of file names assuming they follow the
        Dingo file name format, i.e. Advanced_processed_SITENAME_vXX.csv
        """
        file_name = [ re.compile('^\w+').findall(f) for f in ec_files ]
        split_name = [ re.split('_',f) for f in sum(file_name,[]) ]
        name_ver = [ sublist[-2] for sublist in split_name ]
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
            prefix = [ char[:2] for char in sites ]
            return '_'.join(prefix)
        else:
            return sites.pop()

    def approx_optim_force(self, extrema_df, f_mod, p0=[1,2], \
                           xlabel="SWC_smooth", ylabel="NDVI_grass"):
        """
        This function acts as wrapper to fit some arbitrary univariate model given
        X and Y data series. Some dataframe is passed and the two variables of
        interest are extracted based on the two label values, and then optimised
        on. Returns tabular dataframe giving parameter estimates and their errors.
        """
        mo = _mo.model_optim_extras()
        xobs = extrema_df[xlabel]
        yobs = extrema_df[ylabel]
        sig_res = minimize( mo.min_chi2(f_mod, yobs, xobs), p0 )
        sig_err = mo.get_errors(sig_res, len(yobs))
        site_lab = self.create_label(extrema_df["Site"])
        return pd.DataFrame({'Site':site_lab, 'Value':sig_res['x'], 'Error':sig_err})

    def optimize_all_sites(self, dataset):
        """
        Function does 3 separate optimisations to derive the environmental forcing on
        the pendulum:
            [1] ensemble: or all sites are optimised
            [2] in-sample: individual sites are optimised
            [3] out-of-sample: like ensemble but optimises on all sites minus one
        """
        mo = _mo.model_optim_extras()
        # ensemble optimisation
        all_data = pd.concat(dataset)
        all_res = self.approx_optim_force( all_data, mo.sig_mod )
        all_res["Sampling"] = "ensemble"
        # out-of-sample optimisation (horrible syntax has to be a better functional way)
        sites = self.df_pop_site(all_data["Site"])
        out_res = []
        for i,x in enumerate(sites):
            out_res.append( self.approx_optim_force(all_data.query("Site!=x"), \
                                                    mo.sig_mod ) )
            out_res[i]["Site"] = sites[i]
        out_res2 = pd.concat(out_res)
        out_res2["Sampling"] = "out"
        # in-sample optimisation
        ind_res = pd.concat( map( lambda x: self.approx_optim_force(x, mo.sig_mod), \
                                 dataset ) )
        ind_res["Sampling"] = "in"
        # join all tables of optimisation combinations and return
        return pd.concat([all_res,ind_res,out_res2])

