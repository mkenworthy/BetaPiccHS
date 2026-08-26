import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os, sys
import datetime
runtime = os.path.abspath((sys.argv[0])) + " run at " + datetime.datetime.now().strftime("%c")
import betapic_c as bp
from kepler3 import *
from exorings3 import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from PIL import Image
import paths
matplotlib.rcParams.update({'font.size': 17})

if (len(sys.argv)<2):
    print('Error: you must specify a fractional radius of CPD to be searched for.')
    print('    usage: ./python {} 0.30'.format(sys.argv[0]))
    pass
    #sys.exit() #quit()

figsize_a4 = (8.3, 11.7) #a4
figsize_square_a4 = (8.3, 8.3)
figsize_half_a4 = (8.3 , 8.3/1.618)
figsize_2_1_a4 = (8.3, 8.3/2)
figsize_1_2_a4 =  (8.3/1.618, 8.3*0.99)

b = 0.1 #fraction hill sphere for impact parameter

# #vary between the two fractions of the hill sphere
#f_hill = 0.3 # 0.3 or 0.6
t_mid = 58195. * u.day #first primary transit
vmax_tau = 0.0004
vmax_chi = 2
v_max_upper_tau = 0.0002
v_max_mass = 4e19

# #f_hill = 0.3 #0.3 or 0.6
# t_mid = 59415. * u.day #first primary transit
# vmax_tau =0.03
# vmax_chi = 10
# v_max_upper_tau = 0.03
# v_max_mass = 1.1e21 #None

# #print('Compute the circumplanetary disk using a radius of Hill sphere disk in fraction of Hill radius is {:.3f} at transit midpoint of {}'.format(f_hill, t_mid))

#calculate minimal grain size with the assumption that the dust is concentrated 
#at the area mean weighted planetocentric distance of 0.7 $r_{\text{CPD}}$.
rho = 2.5* u.g/u.cm**3
mean_a_constant = (2*10**5 * bp.L_star.value) / ( rho.to(u.kg/u.m**3).value  * bp.M_star.value**(2/3) * bp.M_c.to(u.earthMass).value**(1/3))
mean_a = mean_a_constant * np.sqrt(0.7) *u.micron # from Grant Kennedy, in Section 4

print('constant', mean_a_constant )
print('minimum grain size at 0.7 r_CPD is {:.1e}'.format(mean_a))


#import data from ASTEP, BRITE and bRing
f_bring = Table.read(paths.data /'binned_BRING.dat', format='ascii.ecsv')
f_astep = Table.read(paths.data /'binned_ASTEP.dat', format='ascii.ecsv')
f_brite = Table.read(paths.data /'binned_BRITE.dat', format='ascii.ecsv')
#plt.plot(f_bring['time'][(58195-44 < f_bring['time']) & (f_bring['time'] < 58195+44)], f_bring['flux'][(58195-44 < f_bring['time']) & (f_bring['time'] < 58195+44)])
#correct for bring systematic offset
mask_bring_min = (58195-40*1.5 < f_bring['time']) & (f_bring['time'] < 58195-40)
mask_bring_max = (58195+40 < f_bring['time']) & (f_bring['time'] < 58195+40*1.5)
f_bring['flux'] -= np.mean(np.concatenate([mask_bring_min, mask_bring_max]))

#bring has negative optical depth values, so to stop it contribution to the upper mass limit, it's values are ignored. 
if t_mid == 58195. * u.day:
    f_bring = f_bring[~((58195-40 < f_bring['time'])  & (f_bring['time'] < 58195+40))]

# remove any anomalously low error bars from the BRING data and set them with 1% errors
m=(f_bring['ferr']<0.002)
f_bring['ferr'][m] = 0.01

# calculate Hill sphere for beta pic b
r_hill = rhill(bp.M_star, bp.M_c, bp.a_c).to(u.au) # radius of hill sphere
print('Hill sphere radius {:5.3f}'.format(r_hill))

angle_i = 30    # inclination of the disk in degrees
angle_phi = 80  # tilt of the disk in degrees

impact_b = b*r_hill # impact parameter of the star
print('Impact parameter distance {:5.3f}'.format(impact_b))

vplanet = vcirc(bp.M_star, bp.M_c, bp.a_c)
print('circular velocity at planet is {:5.3f}'.format(vplanet))

# okay, we need a simple function to convert time to x position w.r.t. the planet
# this should be the full kepler orbit, but we'll go with:
def xp(t, t_mid, vtrans):
    return ((t-t_mid)*vtrans).to(u.au)

from lmfit import Model
gmodel = Model(disk_lc_model, independent_vars=['x','xlower','xupper'])
print('parameter names: {}'.format(gmodel.param_names))
print('independent variables: {}'.format(gmodel.independent_vars))

params = gmodel.make_params(deltadisk=0.0, foutside=1.0, xlower=1.8, xupper=3.4)

# number of points in the sampling grid for the plots
n_i = 25
n_phi = 39

i_range = np.linspace(10, 90, n_i)
phi_range = np.linspace(0, 180, n_phi)

def fitdisk(x_coord, f, ferr, i_range, phi_range, r_disk, gmodel, params):
    deltaf = np.ma.zeros((i_range.size,phi_range.size))
    redchisq = np.ma.zeros((i_range.size,phi_range.size))
    deltaf_err = np.ma.zeros((i_range.size,phi_range.size))
    fout = np.ma.zeros((i_range.size,phi_range.size))
    fout_err = np.ma.zeros((i_range.size,phi_range.size))
    sanity = np.ma.zeros((i_range.size,phi_range.size))
    npoints = np.ma.zeros((i_range.size,phi_range.size))
    success = np.ma.zeros((i_range.size,phi_range.size))

    deltaf.mask = True
    redchisq.mask = True
    deltaf_err.mask = True
    fout.mask = True
    fout_err.mask = True
    sanity.mask = True
    npoints.mask = True
    success.mask = True

    for i, curr_i in enumerate(i_range):
        for j, curr_phi in enumerate(phi_range):

            (xring, yring, radring) = sky_to_ring(x_coord, impact_b, curr_i, curr_phi)
            indisk2 = (radring < r_disk)
            nindisk = np.sum(indisk2)
            npoints[i][j] = nindisk
            if (nindisk > 10):

                # find the edges of the disk region
                xp_lower = np.min(x_coord[indisk2])
                xp_upper = np.max(x_coord[indisk2])

                result = gmodel.fit(f, params, x=x_coord, xlower=xp_lower, xupper=xp_upper, weights=1./(ferr))
                deltaf[i][j] = result.params['deltadisk'].value
                redchisq[i][j] = result.redchi
                deltaf_err[i][j] = result.params['deltadisk'].stderr
                fout[i][j] = result.params['foutside'].value
                fout_err[i][j] = result.params['foutside'].stderr
                success[i][j] = result.success

                sanity[i][j] = curr_i

    # propagate the errors for fout and
    # I = I_0 exo(-tau.cos(i))
    # taylor it...
    # I = I_0 * (1 - tau.cos(i))
    # I = I_0 - I_0 tau cos i)
    # tau cos(i) = (I_0 - I)/I0
    # tau cos(i) = 1 - I/I0
    # tau cos(i) = 1 - (fout-deltaf)/fout
    # tau cos(i) = deltaf/fout
    
    tau = deltaf/fout
    tau_err = np.power(np.power((deltaf_err/deltaf),2.)+np.power((fout_err/fout),2),0.5)*tau

    # returning a dict because I'm not sure how much data I'll be returning
    ans_dict = {'redchisq': redchisq, 'df': deltaf, 'dfe':deltaf_err, 'f0':fout, 'f0e':fout_err, 'tau':tau, 'taue':tau_err, 'npoints':npoints, 'success':success}
    return ans_dict

      
def disk_mass(r_disk, tau, mean_a=0.5*u.micron, mean_rho=2.5*u.g/(u.cm*u.cm*u.cm)):
    'simple mass for a face-on circular optically thin disk'
    # cadged from Mellon's derivation in thesis - p.44, eq. 4.9
    Mdisk = (4*np.pi*mean_a*mean_rho*tau*r_disk*r_disk)/3.
    return Mdisk.to(u.g)



#change every time as parameters change
def return_diskfit(f_hill):
    r_disk = r_hill * f_hill # radius of disk
    print('Disk radius is {:5.3f}'.format(r_disk))

    f=f_brite['flux']+1.0
    x_coord = xp(f_brite['time'] * u.day, t_mid, vplanet )
    ferr = f_brite['ferr']
    (brite_disk) = fitdisk(x_coord, f, ferr, i_range, phi_range, r_disk, gmodel, params)

    f=f_bring['flux']+1.0
    x_coord = xp(f_bring['time'] * u.day, t_mid, vplanet )
    ferr = f_bring['ferr']
    (bring_disk) = fitdisk(x_coord, f, ferr, i_range, phi_range, r_disk, gmodel, params)

    f=f_astep['flux']+1.0
    x_coord = xp(f_astep['time'] * u.day, t_mid, vplanet )
    ferr = f_astep['ferr']
    (astep_disk) = fitdisk(x_coord, f, ferr, i_range, phi_range, r_disk, gmodel, params)
    return bring_disk, astep_disk, brite_disk

def plot_diskfit4(d, d2, dataname, i_range, phi_range):
    'plots disk fit results in 2 panels, requires a dict from the disk fitting routine'
    shrink = 0.8
    fig2, f2_axes = plt.subplots(2, 2, figsize=(8.3*1.4, 11.7*0.7), constrained_layout=False)
    (ax1,ax2, ax3, ax4) = f2_axes.flatten()

    # sets the extent=() limits in imshow() so that the centre of each pixel corresponds to i and phi value
    i_min = i_range.min()
    i_max = i_range.max()
    di = 0.5*((i_max-i_min)/(i_range.size-1))

    phi_min = phi_range.min()
    phi_max = phi_range.max()
    dphi = 0.5*((phi_max-phi_min)/(phi_range.size-1))

    im1 = ax1.imshow(d['tau'],
                     extent=(phi_min-dphi,phi_max+dphi,i_min-di,i_max+di),
                     cmap='viridis',
                     origin='lower',
                     vmin=0,
                     vmax = vmax_tau
                    )
    
    im3 = ax3.imshow(d['tau']/d['taue'],
                     extent=(phi_min-dphi,phi_max+dphi,i_min-di,i_max+di),
                     origin='lower',
                     cmap='viridis',
                     vmin=0,
                     vmax = vmax_chi 
                     )
    im2 = ax2.imshow(d2['tau'],
                     extent=(phi_min-dphi,phi_max+dphi,i_min-di,i_max+di),
                     cmap='viridis',
                     origin='lower',
                     vmin=0,
                     vmax = vmax_tau
                    )
    im4 = ax4.imshow(d['tau']/d['taue'],
                     extent=(phi_min-dphi,phi_max+dphi,i_min-di,i_max+di),
                     origin='lower',
                     cmap='viridis',
                     vmin=0,
                     vmax = vmax_chi 
                     )
    cbar = fig2.colorbar(im1, ax=ax1, shrink=shrink)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.minorticks_on()
    cbar3 = fig2.colorbar(im3, ax=ax3, shrink=shrink)   
    cbar3.ax.minorticks_on()
    cbar2 = fig2.colorbar(im2, ax=ax2, shrink=shrink)
    cbar2.formatter.set_powerlimits((0, 0))
    cbar2.ax.minorticks_on()
    cbar4 = fig2.colorbar(im4, ax=ax4, shrink=shrink)
    cbar4.ax.minorticks_on()
    ax1.set_title(r'$\tau\ \sin\ \theta $', fontsize=18)
    ax2.set_title(r'$\tau\ \sin\ \theta $', fontsize=18)
    ax3.set_title(r'Signal to noise of' +'\n' +r'$\tau\ \sin \ \theta$', fontsize=18)
    ax4.set_title(r'Signal to noise of' +'\n' +r'$\tau\ \sin \ \theta$', fontsize=18)
    
    ax3.text(0.05, 0.0, dataname, ha='left', va='bottom', transform=ax3.transAxes) #, **tyb)

    for a in f2_axes.flatten():
        a.xaxis.set_major_locator(MultipleLocator(30))
        a.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        a.yaxis.set_major_locator(MultipleLocator(30))
        a.yaxis.set_major_formatter(FormatStrFormatter('%d'))

        # For the minor ticks, use no labels; default NullFormatter.
        a.xaxis.set_minor_locator(MultipleLocator(10))
        a.yaxis.set_minor_locator(MultipleLocator(10))

        a.set_xlim(-10,190)
        a.set_ylim(-10,100)


    ax3.set_xlabel(r'Tilt $\phi$ [deg]',
                     )
    ax4.set_xlabel(r'Tilt $\phi$ [deg]',
                     )
        
    ax1.set_ylabel(r'Inclination $\theta$ [deg]',
                     )
    ax3.set_ylabel(r'Inclination $\theta$ [deg]',
                     )

    title = fr"Primary transit at MJD {int(t_mid.value )} "
    fig2.suptitle(title,  y=0.95)
    
    plt.tight_layout()
    plt.draw()
    plotout = (paths.figures / 'taumass_{}_diskfit_{}.pdf'.format(t_mid.value,b,t_mid.value,dataname, ))
    plt.savefig(plotout)


bring_disk_03, astep_disk_03, brite_disk_03 = return_diskfit(f_hill=0.3)
bring_disk_06, astep_disk_06, brite_disk_06 = return_diskfit(f_hill=0.6)


#plot_diskfit4(bring_disk_03, bring_disk_06,   'BRING', i_range, phi_range)
#plot_diskfit4(astep_disk_03, astep_disk_06,   'ASTEP', i_range, phi_range)
plot_diskfit4(brite_disk_03, brite_disk_06,   'BRITE', i_range, phi_range)






