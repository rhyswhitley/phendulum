#!/usr/bin/env python

import datetime, time, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
# load own modules
import spring_dynamics as _sd

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

# Plot functions
class model_plotting(object):

    def __init__(self, fig_path):
        """
        Initialise class for all necessary plotting requirements
        variables are set and stored here for easy access and control
        over all plotting
        """
        self.xlabel = "SWC_smooth"
        self.ylabel = "NDVI_norm"
        self.fpath = fig_path
        self.outcol = "#000000"
        self.col = ['#3CB371','#DC143C','#4169E1']
        self.lab = ['ensemble','in','out']


    def plot_allSite_forcing(self, f_mod, data_list, par_list):
        """
        Creates a PDF whose pages contain the results that describes the
        environmental forcing at each site, as well as the forcing based on
        out-of-site sampling and as an ensemble of all sites.
        """
        with PdfPages(self.fpath+'environ_forcing.pdf') as pdf:
            # plot all site points as a reference to the ensemble and out-of-sampling fits
            # plot the optimised site-specific forcing
            for data_i, par_i in zip(data_list, par_list):
                [ self._plot_data(d) for d in data_list  ]
                self._plot_forcing( data_i, par_i, f_mod, pdf )

    def _plot_data(self, data):
        """ Wrapper on a normal plot for preset data illustration """
        plt.plot( data[self.xlabel], data[self.ylabel], 'x', color=self.outcol, alpha=0.3 )

    def _plot_models(self, xs, fxs, color, label ):
        """ Wrapper on a normal plot for preset model illustration """
        plt.plot( xs, fxs, linestyle='-', lw=5, color=color, label=label, alpha=0.5 )

    def _plot_forcing(self, data, k_var, f_mod, pobj):
        """
        Plots the three different expressions of environmental forcing based on
        the three types of sampling used in the analysis
        """
        # Create vectors for model fits
        xs = np.arange(0,0.5,1e-3)
        site_title = set(data["Site"]).pop()

        plt.plot( data[self.xlabel], data[self.ylabel], 'o', color='black', ms=8 )
        [ self._plot_models( xs, f_mod(ks, xs), color=self.col[i], label=self.lab[i] ) for i, ks in enumerate(np.array(k_var)) ]
        plt.xlabel(r'$\theta_{s 10cm}$', fontsize=18)
        plt.ylabel(r'NDVI$_{grass}$')
        plt.axis([0,0.32,-0.05,0.6])
        plt.legend(loc=2)
        plt.title(site_title)
        pobj.savefig()
        plt.close()

    def create_title(self, string):
        """ Splits the site label into two words based on capital letters"""
        return re.sub(r"(\w)([A-Z])", r"\1 \2", string)

    def plot_data_manipulation(self, data_list):
        """
        Creates a PDF whose pages contain the results that describes the
        environmental forcing at each site, as well as the forcing based on
        out-of-site sampling and as an ensemble of all sites.
        """
        with PdfPages(self.fpath+'soilwater.pdf') as pdf:
            [ self._plot_moving_average(d, pdf) for d in data_list ]

        with PdfPages(self.fpath+'ndvi_timeseries.pdf') as pdf:
            [ self._plot_inflexion_points(d, pdf) for d in data_list ]

    def _plot_moving_average(self, data, pobj):
        """
        Adds a Gaussian filter to the soil water content time-seires to
        elucidate major turning points that align with the NDVI time-series.
        """
        xraw = data["SWC10"]
        xmes = data["SWC_smooth"]
        site_title = set(data["Site"]).pop()

        fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot( xraw, '-', color='pink' )
        ax1.plot( xmes, linestyle='-', color='red', lw=2)
        ax2.plot( xraw-xmes, '.', color='pink' )
        ax1.set_ylabel(r'$\theta_{s}$', fontsize=18)
        ax2.set_ylabel(r'$\sigma$', fontsize=18)
        ax2.set_xlabel('Days since ' + data['DT'][0])
        ax1.set_title(site_title)
        pobj.savefig()
        plt.close()

    def _plot_inflexion_points(self, data, pobj):
        yraw = data["NDVI250X"]
        ymes = data["NDVI_grass"]
        yd1 = data["dy1/dt1"]
        yd2 = data["dy2/dt2"]
        site_title = set(data["Site"]).pop()

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(yraw, '-', color='lightblue', label='Total')
        ax1.plot(ymes, '-', color='blue', lw=2)
        ax2.plot(yd1, color='red', lw=2, label="$dx/dt$" )
        ax2.plot(yd2, color='blue', lw=2, label="$d^2x/dt^2$" )
        ax1.set_ylabel('NDVI signal')
        ax2.set_ylabel(r'$x_{t+1}-x_{t}$', fontsize=16)
        ax2.set_xlabel('Days since ' + data['DT'][0])
        ax2.axis([1,len(ymes),-6e-3,6e-3])
        ax2.legend(loc=1)
        ax1.set_title(site_title)
        pobj.savefig()
        plt.close()

#================================================================================
# Pendulum
#================================================================================

    def plot_allSite_pendulum(self, data_list, par_list, f_mod):
        """
        Creates a PDF whose pages contain the results that describes the
        environmental forcing at each site, as well as the forcing based on
        out-of-site sampling and as an ensemble of all sites.
        """
        with PdfPages(self.fpath+'phendulum.pdf') as pdf:
            # plot all site points as a reference to the ensemble and out-of-sampling fits
            # plot the optimised site-specific forcing
            for data, par in zip(data_list, par_list):
                self._plot_pendulum( data, par, f_mod, pdf )

    def _plot_pendulum(self, data, kvar, f_mod, pobj):

        # now assign the optimised coefficients to the pendulum and calculate the motion
        y_grass = data["NDVI_norm"]
        x_mes = data["SWC10"]

        # based on the number of samplings get the prediction of motion
        springs = [ _sd.spring(k, x_mes, f_mod, x_init=0.2, v_init=0.02) for k in np.array(kvar) ]
        y_mod = [ p.calc_dynamics()['x'] for p in springs ]
        force_mod = [ f.calc_dynamics()['Fe'] for f in springs ]
        accel_mod = [ a.calc_dynamics()['a'] for a in springs ]
        veloc_mod = [ v.calc_dynamics()['v'] for v in springs ]

        site_title = set(data["Site"]).pop()

        fig = plt.figure(figsize=(10,8))
        # setup plotting grid
        gs = gridspec.GridSpec(3, 1, height_ratios=[2,1,1])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        # remove xlabels on the second and third plots
        plt.setp(ax1.get_xticklabels(), visible=False)
        # plot data
        ax1.plot( y_grass, color='black', lw=2, label="MODIS" )
        [ ax1.plot( y_mod[i], lw=2, alpha=0.8, color=self.col[i], label=self.lab[i] ) for i, ks in enumerate(y_mod) ]
        ax2.plot( veloc_mod, color=self.col[2], alpha=0.8, lw=1.5 )
        ax3.plot( accel_mod, color=self.col[1], alpha=0.8, lw=1.5 )
        # labels
        ax1.set_ylabel( r"NDVI", fontsize=14 )
        ax2.set_ylabel( r"$\delta$NDVI/$\delta$t", fontsize=18 )
        ax2.set_ylabel( r"$\delta^{2}$NDVI/$\delta$t$^{2}$", fontsize=18 )
        # legends
        ax1.legend(loc=1)
        ax1.set_title(site_title, size=20)

        gs.tight_layout(fig, rect=[0, 0, 1, 0.97])

        pobj.savefig()
        plt.close()

