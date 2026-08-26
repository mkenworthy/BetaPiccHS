"""Plot the 1981 data and the fit"""
from astropy import constants as c
from astropy import units as u
import numpy.ma as ma
import numpy as np
from matplotlib.collections import PatchCollection
from astropy.table import Table
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io import ascii
import betapic_c as bp
import os, sys
import datetime
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import paths

runtime = os.path.abspath((sys.argv[0])) + " run at " + datetime.datetime.now().strftime("%c")

# choose epochs where to fit the 1981 function
fwhm_lamers = 3.2 #days is estimated from Lamers 1997 A&A 328 321 page 8 Figure 7

def m1981(t, t0, peak, bgnd, fwhm=fwhm_lamers, inner_width=0.25, depth=-0.009):
    """m1981 - a model for the 1981 event
    modelled with two components:
    1. a gaussian function with amplitude of `peak` and FWHM of `fwhm`
    2. narrow triangular absorption trough at the midpoint

    t - sampling points for the function
    t0 - the epoch of the central peak
    peak - amplitude of the central peak
    bgnd - the background flux level
    fwhm - full width half max of the gaussian curve
    inner_width - width of the narrow eclipser
    depth - relative depth of the narrow eclipser"""

    dt = (t-t0)

    # make the gaussian function
    
    # FWHM = 2.sqrt(2 ln 2)sig
    sig = fwhm / 2.355
    di = peak*np.exp(-dt*dt/(2*sig*sig))

    # mask central peak and replace with narrow eclipser
    mask = np.abs(dt)<inner_width
    di_edge = peak*np.exp(-inner_width*inner_width/(2*sig*sig))
    # y = mx + c
    # dt = 0, di = depth
    # dt = inner_width, di = di_edge

    m = (di_edge - depth)/(inner_width)
    di[mask] = depth + m*np.abs(dt[mask])

    di = di + bgnd
    return(di)

# Lecavelier des Etangs photometry
# leclavier des etangs 1992 AA 328 311 - Table 1
# beta pic photometry
t_lde = Table.read( """     JD          Vmag
                            4914.780    3.834
                            4914.857    3.836
                            4917.804    3.824
                            4917.857    3.824
                            4918.628    3.805
  #                          4918.720    3.835
                            4918.786    3.838
                            4918.856    3.845
                            4919.802    3.823
                            4919.853    3.824
                            4920.787    3.828
                            4920.859    3.828
                            4925.791    3.839
                            4925.847    3.839
                    """, format='ascii')

# The complete beta pic photometry from Lecavelier 1995
t = ascii.read(paths.data / 'lecavelierdesetangs1995/table', format='cds', 
               readme=paths.data / 'lecavelierdesetangs1995/ReadMe')
t_1981epoch = t['JD'] - 2440000.

fig, ax1 = plt.subplots(1,1,figsize= (8.3*0.5 , 8.3/1.618*0.6), constrained_layout=True)

# Lecavelier 1995 photometry
ax1.scatter(t_1981epoch, t['Vmag'], color='grey', s=10)

t_mid = 4919.04 # from Lecavelier des Etangs 1997
t_mid = t_mid - 0.14 # seems to be an offset I need by looking at the Lamers 1997 Figure 7

V_sigma          = 0.005 * np.ones_like(t_lde['JD']) # error quoted in Lamers 1997 Figure 1
V_mag_background = 3.842 # V band mean magnitude from Lamers 1997 Figure 1 estimate
V_1981_peak      = 0.034 # Amplitude of the broad peak model from Lamers 1997 estimated from Figure 7

ax1.errorbar(t_lde['JD'], t_lde['Vmag'], yerr=V_sigma,
             fmt='o', color='red',ecolor='red',capsize=0 ,mew=1, elinewidth=1,ms=2)
ax1.set_xlabel('MJD [days]')
ax1.set_ylabel('V band [mag]')

dt = 8. #half width of the figure plot
ax1.set_ylim(3.86,3.78)
ax1.set_xlim(t_mid-dt, t_mid+dt)

t = np.arange(t_mid-dt, t_mid+dt, 0.05)
ax1.xaxis.set_minor_locator(MultipleLocator(1))  
ax1.yaxis.set_minor_locator(MultipleLocator(0.005))   
ax1.plot(t, m1981(t, t_mid, -V_1981_peak, V_mag_background, fwhm=fwhm_lamers, depth=0.009)) 
plt.draw()
plt.savefig(paths.figures / '1981_eventm1981model.pdf', bbox_inches='tight')
