#!/usr/bin/env python2

import datetime, time
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats.mstats import mquantiles
# own modules
import model_optim_extras as _mo
import spring_dynamics as _sd
import data_handling as _dh

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

def main():
    springObj = import_pickle(mcfile)
    _plotPosteriors(springObj)
    _enviroForce(springObj)
    _plotSpring(springObj)

def import_pickle(fpath):
    fileObject = open(fpath, 'rb')
    data = pickle.load(fileObject)
    fileObject.close()
    return data

def _plotPosteriors(pymc_obj):

    k_samples = [ pymc_obj['k_{0}'.format(i)].values() for i in range(5)]
    noSamples = len(k_samples)

    fig = plt.figure(figsize=(8,9))
    gs = gridspec.GridSpec(noSamples, 2, width_ratios=[2,1])
    for i,j in zip(range(0,10,2), range(noSamples)):
        _myTrace(fig, gs[i], k_samples[j][0])
        _myHist(fig, gs[i+1], k_samples[j])
    plt.show()

def _myHist(fObj, gObj, sample):
    pObj = fObj.add_subplot(gObj)
    pObj.hist(sample, histtype='stepfilled', color="#DC143C", edgecolor="none",
            normed=False, bins=35, alpha=0.8)

def _myTrace(fObj, gObj, sample):
    pObj = fObj.add_subplot(gObj)
    pObj.plot(sample, color="#DC143C", linestyle='-', alpha=0.8)

def _enviroForce(pymc_obj):
    k_samples = [ pymc_obj['k_{0}'.format(i)].values() for i in range(5)]
    k_means = [ np.mean(k) for k in k_samples ]

    xs = np.arange(0., 0.3, 0.001)
    eforce = mo.sig_mod2( k_means[0:3], xs )

    plt.plot( xs, eforce, lw=3 )
    plt.show()

def _plotSpring(pymc_obj):
    k_samples = [ pymc_obj['k_{0}'.format(i)].values() for i in range(5)]
    k_means = [ np.mean(k) for k in k_samples ]
    k_quantU = [ np.percentile(k, 97.5) for k in k_samples ]
    k_quantL = [ np.percentile(k, 2.5) for k in k_samples ]


    # import data as a list of Pandas dataframes
    dat_path = "../data/"
    raw_data = dh.import_data(dat_path)
    # map data transformations to each dataset imported
    cor_data = [ dh.grass_correct_data(rd) for rd in raw_data ]

    prediction = _sd.spring(k_means, cor_data[5]["SWC10"], mo.sig_mod2) \
                    .calc_dynamics()['x']
    pred_ci_U = _sd.spring( k_quantU, cor_data[5]["SWC10"], mo.sig_mod2) \
                    .calc_dynamics()['x']
    pred_ci_L = _sd.spring( k_quantL, cor_data[5]["SWC10"], mo.sig_mod2) \
                    .calc_dynamics()['x']

    #plt.fill_between( range(len(prediction)), pred_ci_U,  pred_ci_L, color='red', alpha=0.2)
    plt.plot( cor_data[5]["NDVI_grass"], '.', lw=2, color="black")
    plt.plot( prediction, lw=2, color="#DC143C")
    plt.show()

if __name__=="__main__":
    mcfile = "../outputs/spring_trace"
    mo = _mo.model_optim_extras()
    dh = _dh.data_handling()
    main()

