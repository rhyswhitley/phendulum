#!/usr/bin/env python2.7

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import colorsys, sys
import socket
import calendar
import datetime
from matplotlib.backends.backend_pdf import PdfPages
from os.path import expanduser
#from matplotlib.colors import LinearSegmentedColormap
#import scipy.integrate as integrate

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

class retrieval_system(object):
    """
    This class simply imports a large eddy covariance data set for analysis using the pandas library
    while performing rudimentary error checking.

    get_data() - returns the data to the user given a correct file path
    import_data() - imports the CSV file through pandas
    echo_variables() - prints the column names and their index to the screen
    """

    def __init__(self, file_path, verbose=False):
        # create a ASCII color class
        self.bc = bcolors
        # functor for parse datetime import
        self.dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        # NULL default file path
        self.fpath = file_path
        self.verbose = verbose
        return

    def get_data(self):
        # attempt to import the data
        try:
            # import data
            self.import_data()
            if self.verbose == True:
                print "File located in {}".format(self.fpath)
            # return column names to the user
            self.get_variables()
            # return imported dataset to user
            return self.ecdata
        except:
            print self.bc.FAIL + "*** NO DATA ***"
            return None
        return

    def import_data(self):
        print "Connected to: " + socket.gethostname() + "\n"
        try:
            # read in data
            self.ecdata = pd.read_csv( self.fpath,  parse_dates=True , index_col=['DT'])
            print self.bc.OKGREEN + ">> File loaded.\n" + self.bc.ENDC
        except IOError:
            print self.bc.FAIL + '>> File not there. Check path!\n'
        except MemoryError:
            print self.bc.FAIL + '>> File too large to be loaded.\n'

    def get_variables(self):
        # now print column names and their indicies
        ncol = self.ecdata.shape[1]
        #self.namelist = pd.DataFrame({'name': self.ecdata.columns.values}, index=range(ncol))
        self.namelist = pd.Series(self.ecdata.columns.values, index=range(ncol))
        if self.verbose == True:
            # print screen headers
            print self.bc.OKBLUE + "Row\tName" + self.bc.ENDC
            for i in range(ncol):
                print "{0}\t{1}".format(self.namelist.index[i], self.namelist[i])

    def search_variables(self, regex):
        # prints to the screen the location (if any) of the variable in question
        try:
            search_var = self.namelist[self.namelist.str.contains(regex)]
            if self.verbose == True:
                print search_var
            return search_var
        except:
            print "variable not present or regex incorrect"
            return None


# File paths that change between computers
site = "AdelaideRiver"
version = "_v11b"
fpath = "~/Dropbox/DINGO OUTPUTS/{0}/Advanced_processed_data_{0}{1}.csv".format(site,version)
# ** if you're not on Rhys' computer then replace the above line with the one below
#fpath2 = "~/Dropbox/ADVANCED RESULTS (1)/{0}/Advanced_processed_data_{0}{1}.csv".format(site,version)
file_loc = expanduser(fpath)

rsys = retrieval_system(file_loc)
data = rsys.get_data()
dvar = rsys.search_variables('^surf*')

# Aggregate values for daily means
gdata = data.groupby([data.index.month,data.index.day]).aggregate(np.mean)

# Resample the water balance quantities as daily sums
precip = data['Precip_Con'].resample('D', how='sum')
etrans = data['Fe_Con'].resample('D', how='sum')
# Take cumulative sums
precip_c = np.cumsum(precip)
etrans_c = np.cumsum(etrans)
# Convert W m-2 to mm d-1
conv = 2.501 / 1800
# Calculate water balance deficit
wbal = etrans_c*conv - precip_c

# waveband length for albedo collections
wb = map(str, [ 648, 858, 470, 555, 1240, 1640, 2130 ] )

def plot_surfrefl(dvar, covar, width=16, height=6):
    n_plots = len(dvar)
    n_col = 3
    n_row = n_plots/3

    # MultIndex isn't always set by tag
    #mindex = gdata.index.get_level_values('Month')
    mindex = gdata.index.get_level_values(0)
    r_colors = cm.rainbow(np.linspace(0, 1, len(set(mindex))))
    m_labels = [calendar.month_name[i] for i in range(1,13)]

    #range_bycol = lambda ni, nc, nr: np.reshape(np.reshape( range(ni), (nc,nr) ).T, ni)
    range_byrow = lambda ni, nc, nr: np.reshape( range(ni), (nc,nr) ).T
    surf_ix = range_byrow(n_plots, n_row, n_col)

    yseq_end = n_col*(n_row - 1) + 1
    xseq_sta = n_plots-n_col + 1
    ylab_loc = np.linspace(1, yseq_end, n_row, endpoint=True).tolist()
    xlab_loc = np.linspace( xseq_sta, n_plots, n_col, endpoint=True).tolist()

    r_colors = cm.rainbow(np.linspace(0, 1, len(set(mindex))))
    m_labels = [calendar.month_name[i] for i in range(1,13)]

    fig = plt.figure(figsize=(width, height))
    sax = fig.add_subplot(111) # <<< big subplot
    # Turn off axis lines and ticks of the big subplot
    sax.spines['top'].set_color('none')
    sax.spines['bottom'].set_color('none')
    sax.spines['left'].set_color('none')
    sax.spines['right'].set_color('none')
    sax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    gs = gridspec.GridSpec(n_col, n_row)
    for i in range(n_col):
        for j in range(n_row):
            ax = fig.add_subplot(gs[i, j], axisbg="gray")
            for k, c ,ml in zip(range(12), r_colors, m_labels):
                sdata = gdata[(gdata.index.get_level_values(0)==k+1)]
                pt = plt.scatter(sdata.ix[:,covar], sdata.ix[:,dvar.iloc[surf_ix[i, j]]], \
                        color=c, marker='o', alpha=0.6, label=ml )
            if j==n_row-1 and i==0:
                plt.legend(labels=range(12), loc='upper right', handles=pt, fontsize=8, bbox_to_anchor=(2, 1.05))
            if i+1 in ylab_loc:
                ax.set_xlabel(r'$\lambda_{j}$='+wb[j]+' nm', fontsize=14)
                ax.xaxis.set_label_position('top')
    plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.1, hspace=0.3, wspace=0.4)

plot_surfrefl( dvar, 'Ta_Con' )
sax.set_ylabel('Surface Reflectance', fontsize=14)
sax.set_xlabel(r'Air Temperature ($\degree$C)', fontsize=14)


#with PdfPages('multipage.pdf') as pdf:
#    # Surface Reflectance vs Temperature
#    sax = fig.add_subplot(111) # <<< big subplot
#    # Turn off axis lines and ticks of the big subplot
#    sax.spines['top'].set_color('none')
#    sax.spines['bottom'].set_color('none')
#    sax.spines['left'].set_color('none')
#    sax.spines['right'].set_color('none')
#    sax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
#
#    gs = gridspec.GridSpec(n_col, n_row)
#    for i in range(n_col):
#        for j in range(n_row):
#            ax = fig.add_subplot(gs[i, j], axisbg="gray")
#            for k, c ,ml in zip(range(12), r_colors, m_labels):
#                sdata = gdata[(gdata.index.get_level_values(0)==k+1)]
#                pt = plt.scatter(sdata.ix[:,'Ta_Con'], sdata.ix[:,dvar.iloc[surf_ix[i, j]]], \
#                        color=c, marker='o', alpha=0.6, label=ml )
#            if j==n_row-1 and i==0:
#                plt.legend(labels=range(12), loc='upper right', handles=pt, fontsize=8, bbox_to_anchor=(2, 1.05))
#            if i+1 in ylab_loc:
#                ax.set_xlabel(r'$\lambda$='+wb[j]+' nm', fontsize=14)
#                ax.xaxis.set_label_position('top')
#
#
#    plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.1, hspace=0.3, wspace=0.4)
#    sax.set_xlabel(r'Air Temperature ($\degree$C)', fontsize=14)
#    sax.set_ylabel('Surface Reflectance', fontsize=14)
#    pdf.savefig(facecolor='w')
#    plt.close()
#
#    # Rainfall vs Temperature
#    fig = plt.figure(figsize=(16, 6))
#    sax = fig.add_subplot(111) # <<< big subplot
#    # Turn off axis lines and ticks of the big subplot
#    sax.spines['top'].set_color('none')
#    sax.spines['bottom'].set_color('none')
#    sax.spines['left'].set_color('none')
#    sax.spines['right'].set_color('none')
#    sax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
#
#    gs = gridspec.GridSpec(n_col, n_row)
#    for i in range(n_col):
#        for j in range(n_row):
#            ax = fig.add_subplot(gs[i, j], axisbg="gray")
#            for k, c ,ml in zip(range(12), r_colors, m_labels):
#                sdata = gdata[(gdata.index.get_level_values(0)==k+1)]
#                pt = plt.scatter(sdata.ix[:,'Fc_ustar'], sdata.ix[:,dvar.iloc[surf_ix[i, j]]], \
#                        color=c, marker='o', alpha=0.6, label=ml )
#            if j==n_row-1 and i==0:
#                plt.legend(labels=range(12), loc='upper right', handles=pt, fontsize=8, bbox_to_anchor=(2, 1.05))
#            if i+1 in ylab_loc:
#                ax.set_xlabel(r'$\lambda$='+wb[j]+' nm', fontsize=14)
#                ax.xaxis.set_label_position('top')
#
#
#    plt.subplots_adjust(left=0.05, right=0.9, top=0.9, bottom=0.1, hspace=0.3, wspace=0.4)
#    sax.set_xlabel('Rainfall (mm)', fontsize=14)
#    sax.set_ylabel('Surface Reflectance', fontsize=14)
#    pdf.savefig()
#    plt.close()
#
#    # Set the file's metadata via the PdfPages object:
#    d = pdf.infodict()
#    d['Title'] = 'Surface Reflectance'
#    d['Author'] = u'Rhys Whitley'
#    d['Subject'] = 'Savanna Phenology Analysis'
#    d['Keywords'] = 'savanna grass C4 albedo green'
#    d['CreationDate'] = datetime.datetime(2014, 11, 1)
#    d['ModDate'] = datetime.datetime.today()
#    #sys.exit()
