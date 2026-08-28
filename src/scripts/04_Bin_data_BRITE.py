from astropy import time
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.table import Table, vstack
from matplotlib.colors import LogNorm,ListedColormap
import numpy as np
from astropy.io import fits
import os, sys
import datetime
import paths

from rebinning import *

runtime = os.path.abspath((sys.argv[0])) + " run at " + datetime.datetime.now().strftime("%c")
# add timestamp to plot
tya = dict(color='black', fontsize=8)


#bin BRITE data
fname_brite = paths.data / 'brite/brite_all_R.fits'
data_brite = Table.read(fname_brite)

x_brite = data_brite['time']
y_brite = data_brite['flux']

local_midday = 0.68 # local offset for the middle of the day.
t_plot_start_epoch = np.min(x_brite) - 10.
t_plot_end_epoch = np.max(x_brite) + 10.

t_rebin = np.arange(t_plot_start_epoch, t_plot_end_epoch, 0.05)
(t_brite_binned, f_brite_binned, f_brite_sig_binned, n_poi) = rebintimeseries(x_brite, y_brite, t_rebin)
mb = (f_brite_sig_binned < 0.01)
t_brite_sel, f_brite_sel, f_brite_sig_sel = t_brite_binned[mb], f_brite_binned[mb], f_brite_sig_binned[mb]

#write binned data to file
dat = Table( [t_brite_sel, f_brite_sel, f_brite_sig_sel], names=('time','flux','ferr')  )
dat.write(paths.data / 'binned_BRITE.dat', format='ascii.ecsv', overwrite=True)
