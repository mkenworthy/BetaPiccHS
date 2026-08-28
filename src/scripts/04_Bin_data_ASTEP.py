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


#bin ASTEP data
data_astep = Table.read(paths.data / 'astep_all.fits')
x_astep_all = data_astep['time']
y_astep_all = data_astep['flux']

local_midday = 0.68 # local offset for the middle of the day.
t_plot_start_epoch = np.min(x_astep_all ) - 10.
t_plot_end_epoch = np.max(x_astep_all)+10.


(x_aste, y_aste, day_count_aste) = river2(x_astep_all, y_astep_all, scaler=100)
# select 2017 ASTEP photoemtry mask
(xm2, ym2, mask2) = circsel(x_aste, day_count_aste, yc=57920,yw=80, xc=0.48, xw=0.35)
#ax.scatter(xm2, ym2, color='orange')

# select 2018 ASTEP photoemtry mask
(xm3, ym3, mask3) = circsel(x_aste, day_count_aste, yc=57920+365, yw=80, xc=0.48, xw=0.35)
#ax.scatter(xm3, ym3, color='pink')

x_astep_masked = x_aste[mask2+mask3]
y_astep_masked = y_aste[mask2+mask3]
x_astep_photom_masked = x_astep_all[mask2+mask3]
y_astep_photom_masked = y_astep_all[mask2+mask3]

# epochs for rebinned flux
t_rebin = np.arange(t_plot_start_epoch, t_plot_end_epoch, 0.05)
(t_astep_binned, f_astep_binned, f_astep_sig_binned, n_poi) = rebintimeseries(x_astep_photom_masked, y_astep_photom_masked, t_rebin)

ma = (f_astep_sig_binned < 0.007)
t_astep_sel, f_astep_sel, f_astep_sig_sel = t_astep_binned[ma], f_astep_binned[ma], f_astep_sig_binned[ma]

#write binned data to file
dat = Table( [t_astep_sel, f_astep_sel, f_astep_sig_sel], names=('time','flux','ferr')  )
dat.write(paths.data / 'binned_ASTEP.dat', format='ascii.ecsv', overwrite=True)
