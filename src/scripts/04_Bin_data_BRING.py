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


#Bin bring data
hdul = fits.open(paths.data / 'bring/Reduced_betaPic.fits')
primary = hdul[0].data
cols = hdul[1].columns

data_bRING = hdul[1].data
jd_bRing = data_bRING['jd']
raw_bRing = data_bRING['raw']
eraw_bRing = data_bRING['eraw']
reduced_bRing = data_bRING['reduced']
reducedHF_bRing = data_bRING['reducedHF']

#Get bring data en remove nan
x_bring = jd_bRing-2400000.5 #MJD
m = ~np.isnan(reduced_bRing)
time = x_bring[m]
y_bring_all = np.power(10.,(-reduced_bRing/2.5)) - 1.
y_bring = y_bring_all[m]
x_bring = time

local_midday = 0.68 # local offset for the middle of the day.
t_plot_start_epoch = np.min(x_bring) - 10. #bin for 50 days prior and 50 days after
t_plot_end_epoch = np.max(x_bring) + 10.

# epochs for rebinned flux
t_rebin = np.arange(t_plot_start_epoch, t_plot_end_epoch,0.05) #bin for every 0.05day

# BRING
(t_bring_binned, f_bring_binned, f_bring_sig_binned, n_poi) = rebintimeseries(x_bring, y_bring, t_rebin)
mc = (f_bring_sig_binned < 0.01)
t_bring_sel, f_bring_sel, f_bring_sig_sel = t_bring_binned[mc], f_bring_binned[mc], f_bring_sig_binned[mc]


#write binned data to file
dat = Table( [t_bring_sel, f_bring_sel, f_bring_sig_sel], names=('time','flux','ferr')  )
dat.write(paths.data / 'binned_BRING.dat', format='ascii.ecsv',  overwrite=True)
