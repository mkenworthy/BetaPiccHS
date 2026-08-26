from scipy import stats
import numpy as np
import orbitize
from orbitize import results
from astropy import constants as c
from astropy import units as u
from kepler3 import *
import matplotlib.pyplot as plt
import matplotlib
import astropy
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from astropy.time import Time
import paths

#load orbital bundles

filename=paths.data / "gravity_vandalrv_hgcav2.hdf5"
res = results.Results()
res.load_results(filename)
chain = res.post #posterior output from orbit-fitting process (orbits, varying parameters)
tau_ref_epoch = res.tau_ref_epoch #tau is periastron, closest to star

#search within the range of the Beta pic b campaign
mjds_requested = np.linspace(57700,59500,1100)

# create a new system class from the results. Set some things as dummy variables = 1
sys = orbitize.system.System(2, res.data, 1, 1, mass_err=1, plx_err=1, fit_secondary_mass=True, tau_ref_epoch=res.tau_ref_epoch)
ra, dec, rv = sys.compute_all_orbits(res.post[4], mjds_requested)

#calculate Hill radius pic c
M = 1.75 * u.Msun
Mc = 8.5 * u.Mjup
ac = 2.68 * u.au
ec = 0.208
r_hill_mas = 1000 * (rhill(M, Mc, ac)*(1-ec)).value / 19.44
print("r_hill of Beta pic c is {:.1f} mas".format(r_hill_mas))

#add beta pic b
Mb = 9.00 * u.Mjup
ab = 9.00 * u.au
eb = 0.10

r_hill_b_mas = 1000 * (rhill(M, Mb, ab)*(1-eb)).value / 19.44
print("r_hill b of Beta pic b is {:.1f} mas".format(r_hill_b_mas))

#Create data frame with epoch, mjd of minimum separation and the hill radius

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
fig, ax = plt.subplots(1,1,figsize= (8.3*0.5 , 8.3/1.618*0.6), constrained_layout=True)

epoch_min = [[58000,58490],
             [58500,59000], #ranges of minimum to look into
             [59000,59500]]
df = np.empty([200,3,2]) # pd.DataFrame() #store index, minimum hill radius and its date

for i in np.arange(200):
    ra, dec, rv = sys.compute_all_orbits(res.post[-(i+1)], mjds_requested)
    r = np.sqrt(ra*ra+dec*dec)
    ax.plot(mjds_requested,r[:,2,0]/r_hill_mas,alpha=0.1) #, label = r'\Beta Pic c')
    ax.plot(mjds_requested,r[:,1,0]/r_hill_b_mas,alpha=0.1) #, label = r'Pic b')

    for x, epoch in enumerate(epoch_min):
    
        idx_min = (mjds_requested > epoch[0]) & (mjds_requested< epoch[1])
        min_r = np.min((r[:,2,0])[idx_min])
        df[i,x,0] = min_r/r_hill_mas  #minimum seperation
        df[i,x,1] = mjds_requested[r[:,2,0]==min_r ][0]


#add x-axis for year
time_UTC = lambda mjd: Time(list(mjd), format='mjd').jyear
time_MJD = lambda mjd :Time(list(mjd), format='mjd').mjd
secax = ax.secondary_xaxis('bottom', functions=(time_UTC, time_MJD))
secax.set_xlabel('Epoch [year]')
#plt.grid()

#set minor ticks
secax.xaxis.set_minor_locator(MultipleLocator(0.2))
ax.set_xticks(np.arange(57754, 59580, 365 ))
ax.set_xticks(np.arange(57754-146, 59580, 365/5 ), minor=True)
ax.xaxis.set_ticks_position('top') 
ax.xaxis.set_label_position('top') 

#use the mode of minimum separation
ax.axvline(58194) 
ax.axvline(58707) 
ax.axvline(59414)
ax.axvline(58009, color = 'brown')
#ax.legend(loc = 'upper left')
ax.set_xlabel('Epoch [MJD]')
ax.set_ylabel(r'$\beta$ Pic b and c separation [Hill radius]')

props = dict(boxstyle='round', facecolor='skyblue', alpha=0.5)
ax.text(58194-80, 12.2, 'Primary transit', rotation=90, va='top', fontsize=8, bbox=props) #MJD 58210
ax.text(58707-80, 12.2, 'Secondary transit', rotation=90, va='top', fontsize=8, bbox=props) #MJD 58707
ax.text(59414-80, 12.2, 'Primary transit', rotation=90, va='top', fontsize=8,bbox=props) #MJD 59413
ax.text(58009-90, 12.2,  r'$\beta$ Pic b transit', rotation=90, va='top',  fontsize=8, bbox=props) #MJD 58009
plt.xlim(57650, 59500)
#plt.tight_layout()
plt.savefig(paths.figures / 'orbitize_b_c.pdf')
plt.close()

#print midpoints of transits and their error
for i in range(3):
    mode_mjds= stats.mode(df[:,i,1])[0]
    print(f'For epoch {i} we get {mode_mjds :.2f} +{np.max(df[:,i,1])-mode_mjds :.2f} -{mode_mjds-np.min(df[:,i,1]) :.2f} days')

