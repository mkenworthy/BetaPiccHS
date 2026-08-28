import numpy as np
import orbitize
from orbitize import results
from astropy import constants as c
from astropy import units as u
from kepler3 import *
import matplotlib.pyplot as plt
import matplotlib

import paths

matplotlib.rcParams.update({'font.size': 16})
figsize_a4 = (8.3, 11.7) #a4
figsize_square_a4 = (8.3, 8.3)
figsize_half_a4 = (8.3 , 8.3/1.618)
figsize_2_1_a4 = (8.3, 8.3/2)

filename=paths.data / "withc_gravity_vandalrv_11eps.hdf5"

res = results.Results()
res.load_results(filename)
chain = res.post #posterior output from orbit-fitting process (orbits, varying parameters)
tau_ref_epoch = res.tau_ref_epoch #tau is periastron, closest to star

sample_number = 10

mjds_requested = np.linspace(57800,59500,1000)
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
print("r_hill is {:.1f} mas".format(r_hill_mas))

#add beta pic b

Mb = 9.00 * u.Mjup
ab = 9.00 * u.au
eb = 0.10

r_hill_b_mas = 1000 * (rhill(M, Mb, ab)*(1-eb)).value / 19.44
print("r_hill b is {:.1f} mas".format(r_hill_b_mas))


#Create data frame with epoch, mjd of minimum separation and the hill radius

P = 1121
P_err = 15
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
fig, ax = plt.subplots(1,1,figsize=figsize_half_a4)
print("r_hill b is {:.1f} mas".format(r_hill_b_mas))

epoch_min = [[200,500,900], #ranges of minimum to look into
             [400,700,1100]]
df = pd.DataFrame() #store index, minimum hill radius and its date

for i in np.arange(200):
    ra, dec, rv = sys.compute_all_orbits(res.post[i], mjds_requested)
    r = np.sqrt(ra*ra+dec*dec)
    ax.plot(mjds_requested,r[:,2,0]/r_hill_mas,alpha=0.1) #, label = r'\Beta Pic c')
    
    for x in np.arange(3):
        #store index, minimum hill radius and its date
        min_radius = np.min(r[:,2,0][epoch_min[0][x]:epoch_min[1][x]]/r_hill_mas)
        min_mjds = mjds_requested[np.where(r[:,2,0]/r_hill_mas == min_radius)][0]
        df_new_row =pd.DataFrame({"n":[x],"mjds":[min_mjds],"min_r":[min_radius]}) 
        df = pd.concat([df, df_new_row], ignore_index=True)
    ax.plot(mjds_requested,r[:,1,0]/r_hill_b_mas,alpha=0.1) #, label = r'Pic b')

#use the mode of minimum separation
ax.axvline(58210) 
ax.axvline(58707) 
ax.axvline(59413)
ax.axvline(58009, color = 'brown')
ax.set_xlabel('Epoch [MJD]')
ax.set_ylabel(r'$\beta$ Pic b and c separation [Hill radius]')
ax.text(58210-55, 12.4, 'Primary transit at MJD 58210', rotation=90, va='top', fontsize=13)
ax.text(58707-55, 12.4, 'Secondary transit at MJD 58707', rotation=90, va='top', fontsize=13)
ax.text(59413-55, 12.4, 'Primary transit at MJD 59413', rotation=90, va='top', fontsize=13)
ax.text(58009-55, 12.4,  r'$\beta$ Pic b transit at MJD 58009', rotation=90, va='top', fontsize=13)
plt.savefig(paths.figures / 'orbitize_b_c.pdf')
#plt.show()
