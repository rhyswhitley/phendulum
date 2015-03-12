#!/usr/bin/env python2

import pandas as pd
import datetime, time
import pickle, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats import kde
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
    klen = len([ key for key in springObj.keys() if re.match("^k.*",key) ])
    _plotPosteriors(springObj, klen)
    _plotForce(springObj, klen)
    _plotSpring(springObj, klen)
    #_plotCorrelations(springObj, klen)

def import_pickle(fpath):
    fileObject = open(fpath, 'rb')
    data = pickle.load(fileObject)
    fileObject.close()
    return data

def _plotPosteriors(pymc_obj, klen):

    k_samples = [ pymc_obj['k_{0}'.format(i)].values() for i in range(1,klen+1) ]
    noSamples = len(k_samples)

    fig = plt.figure(figsize=(8,9))
    gs = gridspec.GridSpec(noSamples, 2, width_ratios=[2,1])
    for i, j in zip(range(0,8,2), range(noSamples)):
        _plotTrace(fig, gs[i], k_samples[j][0][::1])
        _plotHist(fig, gs[i+1], k_samples[j][::1])
    plt.show()

def _plotCorrelations(pymc_obj, klen):
    k_samples = [ pymc_obj['k_{0}'.format(i)].values() for i in range(1,klen+1)]
    noSamples = len(k_samples)

    fig = plt.figure(figsize=(11,8))
    gs = gridspec.GridSpec(noSamples, noSamples)
    for i in range(noSamples):
        for j in range(noSamples):
            if i==j:
                ax = fig.add_subplot(gs[i,j])
                ax.plot(range(10), color='k')
            elif i<j:
                ax = fig.add_subplot(gs[i,j])
                ax.scatter(k_samples[i][0], k_samples[j][0], color='k', alpha=0.05)
            else:
                _plotHist2D(fig, gs[i,j], k_samples[i][0], k_samples[j][0])
                #_plotKDE2(fig, gs[i,j], k_samples[i][0], k_samples[j][0])
    #plt.tight_layout()
    plt.show()

def _plotKDE1(fObj, gObj, x, y, nbins=30):
    data = np.column_stack((abs(x), abs(y)))
    k = kde.gaussian_kde(data.T)
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    pObj = fObj.add_subplot(gObj)
    pObj.pcolormesh(xi, yi, zi.reshape(xi.shape))

def _plotKDE2(fObj, gObj, x, y, nbins=30):
    #data = np.column_stack((abs(x), abs(y)))
    data = np.append(x, y, axis=1)
    k = kde.gaussian_kde(data.T)
    x_flat = np.r_[x.min():x.max():nbins*1j]
    y_flat = np.r_[y.min():y.max():nbins*1j]
    xi, yi = np.meshgrid(x_flat, y_flat)
    grid_coords = np.append(x.reshape(-1,1), y.reshape(-1,1), axis=1)
    z = k(grid_coords.T)
    z = z.reshape(nbins,nbins)
    pObj = fObj.add_subplot(gObj)
    pObj.imshow(z ,aspect=x_flat.ptp()/y_flat.ptp())

def _plotHist2D(fObj, gObj, sample1, sample2, bins=20):
    pObj = fObj.add_subplot(gObj)
    pObj.hist2d( sample1, sample2, bins=bins)

def _plotHist(fObj, gObj, sample):
    pObj = fObj.add_subplot(gObj)
    pObj.hist(sample, histtype='stepfilled', color="#DC143C", edgecolor="none",
            normed=False, bins=35, alpha=0.8)

def _plotTrace(fObj, gObj, sample):
    pObj = fObj.add_subplot(gObj)
    pObj.plot(sample, color="#DC143C", linestyle='-', alpha=0.8)

def _plotForce(pymc_obj, klen):
    k_samples = [ pymc_obj['k_{0}'.format(i)].values() for i in range(1,klen+1)]
    k_means = [ np.mean(k) for k in k_samples ]

    xs = np.arange(0., 0.3, 0.001)
    #xs = np.arange(0., 1, 0.001)
    eforce = mo.sig_mod1( k_means[0:2], xs )

    plt.plot( xs, eforce, lw=3 )
    plt.show()

def _plotSpring(pymc_obj, klen):
    k_samples = [ pymc_obj['k_{0}'.format(i)].values() for i in range(1,klen+1)]
    k_mean = [ np.mean(k[::10]) for k in k_samples ]
    k_ci05 = [ np.percentile(k[::10], 97.5) for k in k_samples ]
    k_ci95 = [ np.percentile(k[::10], 2.5) for k in k_samples ]

    if klen<5:
        x_init = cor_data.NDVI_norm[0]
        v_init = cor_data.NDVI_norm[1] - cor_data.NDVI_norm[0]
        k_mean = np.concatenate([k_mean,[x_init,v_init]], axis=0)
        k_ci05 = np.concatenate([k_ci05,[x_init,v_init]], axis=0)
        k_ci95 = np.concatenate([k_ci95,[x_init,v_init]], axis=0)

    prediction = _sd.spring(k_mean, cor_data[xlabel], eforce) \
                    .calc_dynamics()['x']
    pred_ci_U = _sd.spring( k_ci05, cor_data[xlabel], eforce) \
                    .calc_dynamics()['x']
    pred_ci_L = _sd.spring( k_ci95, cor_data[xlabel], eforce) \
                    .calc_dynamics()['x']

    plt.fill_between( range(len(prediction)), pred_ci_U,  pred_ci_L, color='red', alpha=0.2)
    plt.plot( cor_data["NDVI_norm"], '.', lw=2, color="black")
    plt.plot( prediction, lw=2, color="#DC143C")
    plt.show()

if __name__=="__main__":
    dh = _dh.data_handling()
    mo = _mo.model_optim_extras()
    eforce = mo.sig_mod1

    xlabel = "SWC10"
    site_name = "AdelaideRiver"
    dat_path = "../data/filtered_"+site_name+"_v12.csv"

    mcfile = "../outputs/mcmc/spring_trace_exp_"+site_name

    raw_data = pd.read_csv(dat_path)
    cor_data = dh.grass_correct_data(raw_data)

    main()

