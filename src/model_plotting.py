#!/usr/bin/env python

import matplotlib.pyplot as plt

# Plot functions
class model_ploting(object):

    def __init__(self, fig_path, dataframe, save_plot=False):
        self.start_date = dataframe.index[0]
        self.xraw = dataframe["SWC10"]
        self.yraw = dataframe["NDVI250X"]
        self.xmes = dataframe["SWC_smooth"]
        self.ymes = dataframe["NDVI_grass"]
        self.yd1 = dataframe["dy1/dt1"]
        self.yd2 = dataframe["dy2/dt2"]
        self.fpath = fig_path
        self.save_plot = save_plot

    def build_report(self):
        """
        Return graphical representations of how environmental forcing has
        been derived.
        """
        self.plot_smooth_swc()

    def is_plotted(self, fig, fname):
        """
        Allows user to switch figure drawing to screen or files.
        """
        if self.save_plot:
            plt.savefig(self.fpath+fname)
            plt.close(fig)
        else:
            plt.show()

    def plot_smooth_swc(self, file_name = "swc_smoothed.pdf"):
        """
        Adds a Gaussian filter to the soil water content time-seires to
        elucidate major turning points that align with the NDVI time-series.
        """
        fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot( self.xraw, '-', color='pink' )
        ax1.plot( self.xmes, linestyle='-', color='red', lw=2)
        ax2.plot( self.xraw-self.xmes, '.', color='pink' )
        ax1.set_ylabel(r'$\theta_{s}$', fontsize=18)
        ax2.set_ylabel(r'$\sigma$', fontsize=18)
        ax2.set_xlabel('Days since 1-Jan-2008')
        self.is_plotted(fig, file_name)

    def plot_grass_predict(self, file_name = "ndvi_corrected.pdf"):
        fig = plt.figure()
        plt.plot( self.yraw, '-', color='lightblue', label='Total')
        plt.plot( self.ymes, '-', color='blue', lw=2, label='Grass' )
        plt.legend(loc=1)
        plt.xlabel('Days since ' + self.start_date)
        plt.ylabel('NDVI')
        self.is_plotted(fig, file_name)

    def plot_inflexion_points(self, file_name="swc_turningpoints.pdf"):
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot( self.ymes, color='purple', lw=2)
        ax2.plot( self.yd1, color='red', lw=2, label="$dx/dt$" )
        ax2.plot( self.yd2, color='blue', lw=2, label="$d^2x/dt^2$" )
        ax1.set_ylabel('NDVI signal')
        ax2.set_ylabel(r'$x_{t+1}-x_{t}$', fontsize=16)
        ax2.set_xlabel('Time-series (t)')
        ax2.axis([1,len(self.ymes),-6e-3,6e-3])
        ax2.legend(loc=1)
        self.is_plotted(fig, file_name)

#================================================================================
# Plot results
#================================================================================

# now assign the optimised coefficients to the pendulum and calculate the motion
model_fitted = sd.spring( xs, txf, spring_opt[0] )
phenology = model_fitted.calc_dynamics()

# plot the results
y_mod = phenology['x']
# there's got to be an easier way than having to convert to arrays
y_res = np.array(y_mod)-np.array(y_grass)

def plot_pendulum():
    fig = plt.figure(figsize=(10,7))
    # setup plotting grid
    gs = gridspec.GridSpec(3, 1, height_ratios=[2.7,1,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax2)
    # remove xlabels on the second and third plots
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    # plot data
    ax1.plot( xs.index, y_grass, color='black', lw=2, label="MODIS" )
    ax1.plot( xs.index, y_mod, color='red', lw=2, label="Pendulum" )
    ax2.plot( xs.index, y_res, '.', color='black', alpha=0.3 )
    ax3.plot( xs.index, xs, color='blue', lw=1.5 )
    # zero line
    ax2.axhline( 0, color='black', linestyle='--' )
    # set axis limits
    ax2.axis([xs.index[0], xs.index[-1], -0.4, 0.4])
    # labels
    ax1.set_ylabel( r"NDVI", fontsize=14 )
    ax2.set_ylabel( r"$\sigma_{NDVI}$", fontsize=18 )
    ax3.set_ylabel( r"$\theta_{s10cm}$", fontsize=18 )
    # legends
    ax1.legend(loc=1)
    ax1.set_title(site, size=20)
    gs.tight_layout(fig, rect=[0, 0, 1, 0.97])
    fig.savefig("../figs/"+site+"_pendulum_fit.pdf")

