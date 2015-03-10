#!/usr/bin/env python

from os import listdir
import re, pickle
import numpy as np
import pandas as pd
# own modules
import data_handling as _dh

def main():
    post_sols = load_data("../outputs/", "^spring_trace_exp_*")

    site_names = dh.get_all_site_names(post_sols[0], regex='\w+$', pos=1)

    # raw traces
    k_traces = [ [ obj['k_{0}'.format(i)].values() for i in range(1,5)] for obj in post_sols[1] ]
    k_means = [ [ summary_stats(t,i,site) for i,t in enumerate(trace) ] for trace,site in zip(k_traces,site_names) ]
    post_table = pd.DataFrame(np.vstack(k_means))
    post_table.columns = ['Site','k','Mean','SD','CI05','CI95']
    post_table["Sampling"] = "in"
    post_table.to_csv("../outputs/spring_posteriors.csv", index=False, index_label=False)

def summary_stats(data,b,site):
    mu = np.mean(data)
    sd = np.std(data)
    ci05 = np.percentile(data, 2.75)
    ci95 = np.percentile(data, 97.5)
    return np.array([site, b, mu, sd, ci05, ci95]).T

def load_data(dat_path, searching):
    """Pass file location and return a list of dictionaries"""
    # get file names
    file_names = [ dat_path+f for f in listdir(dat_path) if re.match(searching,f) ]
    # load pickle object
    pickle_jar = [ read_pickle(fp) for fp in file_names ]
    return (file_names, pickle_jar)

def read_pickle(fpath):
    """Wrapper to read in pickled binary object"""
    fileObject = open(fpath, 'rb')
    data = pickle.load(fileObject)
    fileObject.close()
    return data

if __name__=="__main__":
    dh = _dh.data_handling()
    main()
