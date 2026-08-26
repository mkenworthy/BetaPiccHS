import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import text
from astropy.io import ascii
from astropy.table import Table, vstack
import paths
# BRITE residual light curves are in millimagnitudes with the mean magnitude subtracted.
# The first column is HJD (mid exposure) - 2456000.0 in days.
fname_out = paths.data / 'brite/brite_all_R.fits'

#only use data from BHr, as BTR data has been found to be inconsistent
t1 = ascii.read(paths.data / 'brite/betaPic_2015-16-17-18-BHr.dat')
t2 = ascii.read(paths.data / 'brite/betaPic_2019_BHr.dat')
t3 = ascii.read(paths.data / 'brite/betaPic_2021-BHr-all.dat')
t_b = vstack([t1,t2,t3])

# convert to MJD
time = t_b['col1'] + 2456000.0 - 2400000.5

# the second column is in millimagnitudes
# convert to normalised flux - 1
f1 = np.power(10.,( ((t_b['col2']/1000.)/-2.5))) - 1.

t =data = Table([time, f1], names=('time','flux'))
t.write(fname_out, overwrite=True)
