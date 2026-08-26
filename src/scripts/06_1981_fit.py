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

# # The complete beta pic photometry from Lecavelier 1995
# t = ascii.read('Data/lecavelierdesetangs1995/table', format='cds', 
#                readme='Data/lecavelierdesetangs1995/ReadMe')
# t_1981epoch = t['JD'] - 2440000.

# fig, ax1 = plt.subplots(1,1,figsize= (8.3*0.5 , 8.3/1.618*0.6), constrained_layout=True)

# # Lecavelier 1995 photometry
# ax1.scatter(t_1981epoch, t['Vmag'], color='grey', s=10)

# t_mid = 4919.04 # from Lecavelier des Etangs 1997
# t_mid = t_mid - 0.14 # seems to be an offset I need by looking at the Lamers 1997 Figure 7

# V_sigma          = 0.005 * np.ones_like(t_lde['JD']) # error quoted in Lamers 1997 Figure 1
# V_mag_background = 3.842 # V band mean magnitude from Lamers 1997 Figure 1 estimate
# V_1981_peak      = 0.034 # Amplitude of the broad peak model from Lamers 1997 estimated from Figure 7

# ax1.errorbar(t_lde['JD'], t_lde['Vmag'], yerr=V_sigma,
#              fmt='o', color='red',ecolor='red',capsize=0 ,mew=1, elinewidth=1,ms=2)
# ax1.set_xlabel('MJD [days]')
# ax1.set_ylabel('V band [mag]')

# dt = 8. #half width of the figure plot
# ax1.set_ylim(3.86,3.78)
# ax1.set_xlim(t_mid-dt, t_mid+dt)

# t = np.arange(t_mid-dt, t_mid+dt, 0.05)
# ax1.xaxis.set_minor_locator(MultipleLocator(1))  
# ax1.yaxis.set_minor_locator(MultipleLocator(0.005))   
# ax1.plot(t, m1981(t, t_mid, -V_1981_peak, V_mag_background, fwhm=fwhm_lamers, depth=0.009)) 
# plt.draw()
# plt.savefig('Figures/1981_eventm1981model.pdf', bbox_inches='tight')


"""Make artifical time series and test that it works"""
import numpy as np
import numpy.ma as ma

step = 1.0

def fit_1981(t, f, ferr, t_test_epochs, t_window=8.0, min_npoints=15):
    # t_window - half width of fitting window
    # min_npoints - minimum number of photometric points for a fit within the t_window
    t_test_ampl = np.zeros_like(t_test_epochs) - 1000. # set -1000 to mark bad/missing points
    t_test_ampl_err = np.zeros_like(t_test_epochs) - 1000.

    #
    t_test_b = np.zeros_like(t_test_epochs) - 1000. # set -1000 to mark bad/missing points
    #

    for (i, t_now) in enumerate(t_test_epochs):
        # select the points plus/minus the epoch
    #    print('Trying {:.2f} ...'.format(t_now))
        n_obs_mask = (t>(t_now-t_window)) * (t<(t_now+t_window))
        n = np.count_nonzero(n_obs_mask)
    #    print('{:d} points found'.format(n))

        if n < min_npoints:
            continue

    #    print('nonzero number of points found!')
        t_sel = t[n_obs_mask]
        d_sel = f[n_obs_mask]
        e_sel = ferr[n_obs_mask]

    # add hints and limits to the fit so it doesn't run away
        params = gmodel.make_params(t0=t_now, peak=0.1, bgnd=0.00)
        gmodel.set_param_hint('t0', value=t_now, min=t_now-(step/2.), max=t_now+(step/2.))
        gmodel.set_param_hint('peak', value=0.1, min=0.0, max=5.)
        result = gmodel.fit(d_sel, t=t_sel, bgnd=0.0, t0=t_now, peak=0.1)
        
        if result.success:
            if result.errorbars:
                asdf = result.eval_uncertainty(sigma=3)
                t_test_ampl[i] = result.best_values['peak']
                t_test_ampl_err[i] = result.params['peak'].stderr
                t_test_b[i] = result.best_values['bgnd']
        else:
            print('FAILED to fit at {}'.format(t_now))

        # convert all to masked arrays
        ama = ma.masked_less(t_test_ampl, -999)
        ema = ma.masked_less(t_test_ampl_err, -999)
        tma = np.ma.masked_where(np.ma.getmask(ama), t_test_epochs)
        
        bma = ma.masked_less(t_test_b, -999)
        
    return (tma, ama, ema, bma)


from lmfit import Model
gmodel = Model(m1981, param_names=('t0','peak', 'bgnd'))
print('parameter names: {}'.format(gmodel.param_names))
print('independent variables: {}'.format(gmodel.independent_vars))
params = gmodel.make_params(t0=1050, peak=0.1, bgnd=0)


mlim = 0.035



f_bring = Table.read(paths.data / 'binned_BRING.dat', format='ascii.ecsv')
t_in = np.arange(min(f_bring['time']), max(f_bring['time']), 0.25)

#Bring Data second 
(tmabring, amabring, emabring, bmabring) = fit_1981(f_bring['time'], f_bring['flux'], f_bring['ferr'], t_in)
max_err = 0.05 # too big error bars should be zeroed out
m = (emabring>max_err)
tmabring[m] = ma.masked
amabring[m] = ma.masked
emabring[m] = ma.masked



#Brite data
f_brite = Table.read(paths.data / 'binned_BRITE.dat', format='ascii.ecsv')
(tmabrite, amabrite, emabrite, bmabrite) = fit_1981(f_brite['time'], f_brite['flux'], f_brite['ferr'], t_in)
m = (emabrite>max_err)
tmabrite[m] = ma.masked
amabrite[m] = ma.masked
emabrite[m] = ma.masked

#Astep
f_astep = Table.read(paths.data / 'binned_ASTEP.dat', format='ascii.ecsv')
(tmaastep, amaastep, emaastep, bmaastep) = fit_1981(f_astep['time'], f_astep['flux'], f_astep['ferr'], t_in)
m = (emaastep>max_err)
tmaastep[m] = ma.masked
amaastep[m] = ma.masked
emaastep[m] = ma.masked

tstack = ma.vstack([tmabrite,tmaastep,tmabring])
astack = ma.vstack([amabrite,amaastep,amabring])
estack = ma.vstack([emabrite,emaastep,emabring])
dat = Table( [tstack, astack,estack ], names=('tstack','astack', 'estack')  )
dat.write(paths.data / '1981_fits_Brite_Bring_astep.dat', format='ascii.ecsv', overwrite=True)

min_amp_ind = ma.argmin(estack, axis=0)
min_amp = astack[min_amp_ind,np.arange(min_amp_ind.size)]
min_tim = tstack[min_amp_ind,np.arange(min_amp_ind.size)]
min_err = estack[min_amp_ind,np.arange(min_amp_ind.size)]



fig6, ax = plt.subplots(15, 1, figsize=(8.3, 11.7*1.3), sharex=False, sharey=False, layout="constrained")
ax[0].set_title("Fitting of ASTEP, bRing and BRITE data to the 1981 event") #, fontsize = fs+5)
start = 57770-60
for i in range(15):
    ax[i].set_xticks(np.arange(start+10, start+140,20), np.arange(start+10, start+140,20))
    
    ax[i].set_yticks([0.01, 0.03, 0.05], [0.01, 0.03, 0.05])
    if i ==0:
        ax[i].errorbar(tstack[0], astack[0], yerr=np.abs(estack[0]), fmt='None', color='dodgerblue', 
                       alpha=1, label='BRITE',zorder =10, elinewidth=1)
        ax[i].errorbar(tstack[1], astack[1], yerr=np.abs(estack[1]), fmt='None', color='red', alpha=1, label='ASTEP', elinewidth=1)
        ax[i].errorbar(tstack[2], astack[2], yerr=np.abs(estack[2]), fmt='None', color='green', alpha=1, label='BRING', elinewidth=1)
    else:
        ax[i].errorbar(tstack[0], astack[0], yerr=np.abs(estack[0]), fmt='None', color='dodgerblue', alpha=1, elinewidth=1, zorder =10)
        ax[i].errorbar(tstack[1], astack[1], yerr=np.abs(estack[1]), fmt='None', color='red', alpha=1,  elinewidth=1)
        ax[i].errorbar(tstack[2], astack[2], yerr=np.abs(estack[2]), fmt='None', color='green', alpha=1, elinewidth=1)

        # For the minor ticks, use no labels; default NullFormatter.
    ax[i].xaxis.set_minor_locator(MultipleLocator(5))
    ax[i].tick_params(axis='both', which='major') 
    ax[i].tick_params(axis='both', which='minor') 
    ax[i].set_ylim(-0.006,0.07)
    ax[i].axhline(0, color='black', alpha=0.3)
    ax[i].axhline(mlim, color='black', alpha=0.9, linestyle='dotted')
    ax[i].set_xlim(start, start+150)
    ax[i].set_ylabel("m")
    start +=150
ax[0].legend(loc='center left')

transits = np.array([58195,58707, 59415 ])
transits_min = np.array([20, 7,8 ])
transits_max= np.array([11, 3,8 ])

#color in the transits/eclipse
ax[3].axvspan(58195-22, 58195+22, alpha=0.2 ,color= 'lightskyblue')
ax[3].errorbar(58195, mlim,xerr=[[20],[11]], fmt='o', color="navy" , capsize=10, elinewidth=2, mew =2)

ax[6].axvspan(58707-22, 58707+22, alpha=0.2, color= 'lightskyblue')
ax[6].errorbar(58707, mlim,xerr=[[7],[3]], fmt='o', color="navy", capsize=10, elinewidth=2, mew =2)

ax[11].axvspan(59415-22, 59415+22, alpha=0.2, color= 'lightskyblue')
ax[11].errorbar(59415 , mlim,xerr=[[8],[8]], fmt='o', color="navy", capsize=10, elinewidth=2, mew =2)

#add labels
ax[3].text(0.1,0.8,f'Primary transit at MJD {transits[0]}',   transform=ax[3].transAxes) 
ax[6].text(0.505,0.8,f'Secondary transit at MJD {transits[1]}',   transform=ax[6].transAxes) 
ax[11].text(0.239,0.8,f'Primary transit at MJD {transits[2]}',   transform=ax[11].transAxes) 

ax[3].set_xlabel('MJD [days]') 
ax[14].set_xlabel('MJD [days]') 

    
plt.savefig(paths.figures / '1981_fits.pdf', bbox_inches='tight', dpi=700)
plt.show()
