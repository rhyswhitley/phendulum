#!/usr/bin/env python

from os import listdir
import re, pickle
import numpy as np
import pandas as pd

def main():
    post_sols = load_data("../outputs/", "^spring_trace_exp_*")

    # raw traces
    k_traces = [ [ obj['k_{0}'.format(i)].values() for i in range(1,5)] for obj in post_sols[1] ]
    k_means = [ [ summary_stats(t) for t in trace ] for trace in k_traces ]
#    k_stand = [ [ np.std(t) for t in trace ] for trace in k_traces ]
#    k_ci05 = [ [ np.percentile(t, 2.75) for t in trace ] for trace in k_traces ]
#    k_ci95 = [ [ np.percentile(t, 97.5) for t in trace ] for trace in k_traces ]
    print pd.DataFrame(np.vstack(k_means))

def summary_stats(data):
    mu = np.mean(data)
    sd = np.std(data)
    ci05 = np.percentile(data, 2.75)
    ci95 = np.percentile(data, 97.5)
    return (mu, sd, ci05, ci95)

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
    main()
