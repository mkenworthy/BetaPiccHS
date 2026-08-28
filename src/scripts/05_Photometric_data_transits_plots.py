from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib as mpl
import paths
#load binned data
f_bring = Table.read(paths.data /'binned_BRING.dat', format='ascii.ecsv')
f_astep = Table.read(paths.data / 'binned_ASTEP.dat', format='ascii.ecsv')
f_brite = Table.read(paths.data / 'binned_BRITE.dat', format='ascii.ecsv')

transit_times = np.array([58195,58707, 59415 ])
epoch1 = 58195 #+11.46 -19.65 days
epoch2 = 58707 #+3.28 -6.55 days
epoch3 = 59415 #+8.19 -8.19 days
color = ['r', 'green', 'dodgerblue']
name_instrument = ["BRING", "ASTEP", "BRITE"]

#err_epoch = 100
d_lim =80 #limit transit end and start in days
data_sets = [f_bring, f_astep, f_brite ]

fig = plt.figure(figsize=(8.3 , 8.3/1.618*1.1)) 
subfigs = fig.subfigures(nrows=3, ncols=1, hspace=0.1)
axs = [subfig.subplots(nrows=3, ncols=1, gridspec_kw={'hspace': 0.01}) for subfig in subfigs.ravel()]
for x, subax in enumerate(axs):
    for y, ax in enumerate(subax):
        data = data_sets[y]
        ax.errorbar(data['time'], data['flux'], yerr=data['ferr'], fmt='.', alpha=0.5, label =name_instrument[y], 
                    color =color[y],capsize=1, markersize=1)
        ax.text(0.92,0.1,name_instrument[y], transform=ax.transAxes,  fontweight='bold') 
        ax.errorbar(epoch1, 0,xerr=[[20],[11]], fmt='o', color="navy", capsize=10, elinewidth=1, mew =1)
        ax.errorbar(epoch2, 0,xerr=[[7],[3]], fmt='o', color="navy", capsize=10, elinewidth=1, mew =1)
        ax.errorbar(epoch3, 0,xerr=[[8],[8]], fmt='o', color="navy", capsize=10, elinewidth=1, mew =1)
        ax.fill_between(np.linspace(epoch1-22, epoch1+22) , -0.2, 0.2,
                     color = 'lightskyblue', alpha =0.2)
        ax.fill_between(np.linspace(epoch2-22, epoch2+22) , -0.2, 0.2,
                     color = 'lightskyblue', alpha =0.2)
        ax.fill_between(np.linspace(epoch3-22, epoch3+22) , -0.2, 0.2,
                     color = 'lightskyblue', alpha =0.2)
        ax.axvline(transit_times[x], color='lightskyblue', linestyle = '--')
        
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        if y ==0:
            ax.set_ylim(-0.059,0.059)
            ax.yaxis.set_major_locator(MultipleLocator(0.035))
            ax.set_xticks([])
        if y ==1:
            ax.set_ylim(-0.029,0.029)
            ax.yaxis.set_major_locator(MultipleLocator(0.02))
            ax.set_xticks([])
            ax.set_ylabel("Normalized flux")
        if y ==2:
            ax.set_ylim(-0.013,0.013)
            ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        if x==0 and y==0:
            ax.set_title(f"Primary transit at MJD {epoch1}" )
        if x==1 and y==0:
            ax.set_title(f"Secondary transit at MJD {epoch2}")
        if x==2 and y==0:
            ax.set_title(f"Primary transit at MJD {epoch3}")
        ax.set_xlim(transit_times[x]-d_lim , transit_times[x]+d_lim)
        if x==2 and y==2:
            ax.set_xlabel('Epoch [MJD ]')

plt.savefig(paths.figures / 'All_photometric_transits.pdf', bbox_inches='tight')
