#!/usr/bin/env python

import numpy as np
import pandas as pd
from os import listdir
import re

def main(fpath, site):

    # echo to user
    print("Creating table for site => "+ site)

    # Display column names and their index for reference
    ec_data = pd.read_csv(fpath, parse_dates=True , index_col=['DT'])

    # Let just grab what we need -- SWC and NDVI
    _map = ec_data.groupby(ec_data.index.year)["Precip_Con"].sum()
    _mat = ec_data.groupby(ec_data.index.year)["Ta_Con"].mean()

    return (site, _map.mean(), _mat.mean(), _map.std(), _mat.std())


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

def whittle_files( _fpath, _names ):
    # get only csv files
    get_files = [ f for f in listdir(_fpath) if f.endswith('.csv') ]
    # take only what is in the list of names
    grab_these = pd.Series(get_files).str.contains("|".join(_names))
    return pd.Series(get_files)[np.array(grab_these)]


if __name__ == '__main__':

    # Import EC tower Dingo dataset for Sturt Plains [move this to a config file]
    out_path = "../data/"
    fig_path = "../figs/"
    search_path = out_path+"Dingo_v12/"
    site_path = search_path+"site coordinates ozflux.csv"
    show_names = False
    draw_plot = True
    out_name = "site_char.csv"
    version = "_v12"

    geocord = pd.read_csv(site_path)

    # collect all files for processed eddy covariance datasets in the data folder
    natt_names = ["AdelaideRiver","AliceSprings","DalyUncleared","DryRiver", \
                  "HowardSprings","SturtPlains"]

    natt_files = whittle_files(search_path, natt_names)
    natt_path = [ search_path+f for f in natt_files ]
    site_name = [ get_site_name(f) for f in natt_files ]
    res = [ main(p,s) for p,s in zip(natt_path,site_name) ]
    clim = pd.DataFrame(res)
    clim.columns = ["Label","MAP","MAT","sig_MAP","sig_MAT"]

    table_out = pd.merge(clim, geocord, left_on='Label', right_on='Label' )

    table_out.to_csv(out_path+out_name, sep=",", index=False, index_label=False)
