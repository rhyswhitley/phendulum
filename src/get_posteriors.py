#!/usr/bin/env python

from os import listdir
import re, pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
# own modules
import data_handling as _dh

def plot_posteriors(traces, burn, lag, site_names, blen=20):

    n_sites = len(traces)
    plab = ["slope","bias","resist","drag"]

    fig = plt.figure(figsize=(9,9))
    gs = gridspec.GridSpec(n_sites, 4)

    sax = fig.add_subplot(111) # <<< big subplot
    # Turn off axis lines and ticks of the big subplot
    sax.spines['top'].set_color('none')
    sax.spines['bottom'].set_color('none')
    sax.spines['left'].set_color('none')
    sax.spines['right'].set_color('none')
    sax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    for i in range(n_sites):
        for j in range(4):

            pObj = fig.add_subplot(gs[i,j])

            pObj.tick_params(axis='y', which='both', labelleft='off')
            if i<n_sites-1:
                pObj.tick_params(axis='x', which='both', labelbottom='off')

            trace = traces[i][j][0]
            sample = trace[burn:len(trace):lag]
            p_true = sample.mean()
            p_ci05 = np.percentile(sample, 2.75)
            p_ci95 = np.percentile(sample, 97.5)
            binny = np.histogram(sample, bins=blen, normed=True, density=False)
            hmax = max(binny[0])*1.1
            pObj.hist( sample, histtype='stepfilled', color="#DC143C", \
                    edgecolor="black", normed=False, bins=blen, alpha=0.2)
            pObj.vlines(p_true, 0, hmax, linestyle="-", lw=2, label="$\mu$")
            pObj.vlines(p_ci05, 0, hmax, linestyle="--", lw=1.5, label="CI 05%", alpha=0.5)
            pObj.vlines(p_ci95, 0, hmax, linestyle="--", lw=1.5, label="CI 95%", alpha=0.5)
            handles, labels = pObj.get_legend_handles_labels()

            if i==n_sites-1:
                pObj.set_xlabel('$k_{'+plab[j]+'}$', size=16)
            if j==0:
                pObj.set_ylabel(site_names[i], size=12)
                pObj.axis([0,100,0,hmax])
            if j==1:
                pObj.axis([0,8,0,hmax])
            if j==2:
                pObj.axis([0,2,0,hmax])
            if j==3:
                pObj.axis([0,5,0,hmax])

    sax.legend(handles, labels, loc='center', ncol=3, bbox_to_anchor=(0.5, -0.1), prop={'size':11} )
    plt.tight_layout(h_pad=0.5, w_pad=0)
    plt.savefig("../figs/posteriors_distributions.pdf")

def plot_traces(traces, lag):
    n_sites = len(traces)
    r_colors = cm.rainbow(np.linspace(0, 1, n_sites))

    fig = plt.figure(figsize=(8,9))
    gs = gridspec.GridSpec(n_sites,4)

    for i in range(n_sites):
        for j in range(4):
            pObj = fig.add_subplot(gs[i,j])
            sample = traces[i][j][0]
            pObj.plot(sample[::lag], color=r_colors[i], alpha=0.8)
    plt.tight_layout(h_pad=1, w_pad=1)
    plt.show()

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

def main():
    post_sols = load_data("../outputs/mcmc/", search)

    site_names = dh.get_all_site_names(post_sols[0], regex='\w+$', pos=1)

    # raw traces
    k_traces = [ [ obj['k_{0}'.format(i)].values() for i in range(1,5)] \
                for obj in post_sols[1] ]
    k_means = [ [ summary_stats(t,i,site) for i,t in enumerate(trace) ] \
                for trace,site in zip(k_traces,site_names) ]
    post_table = pd.DataFrame(np.vstack(k_means))
    post_table.columns = ['Site','k','Mean','SD','CI05','CI95']
    post_table["Sampling"] = "in"
    #post_table.to_csv("../outputs/pars/spring_posteriors.csv", index=False, index_label=False)
    plot_posteriors(k_traces, burn, lag, site_names)
    #plot_traces(k_traces, lag)

if __name__=="__main__":
    dh = _dh.data_handling()
    lag = 5
    burn = 5000
    #search = "^spring_trace_uniInit_*"
    search = "^spring_trace_exp_*"
    main()
