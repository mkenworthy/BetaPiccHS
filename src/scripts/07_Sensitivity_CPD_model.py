import numpy as np
import matplotlib.pyplot as plt
import betapic_c as bp
from kepler3 import *
from exorings3 import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from lmfit import Model
import paths
#vary this parameter for the mock data to get the CPD of different offsets from the set midpoint of the transit
offset_days = np.array([0,3,5,10,15,20])

import os, sys
import datetime
runtime = os.path.abspath((sys.argv[0])) + " run at " + datetime.datetime.now().strftime("%c")

# for a given tip and tilt, you can calculate the times where you cross that chord across the disk
# if you do hit that radius, you then split the photometry into 'yes we hit it' and 'no we don't'.

# then fit for the flux change between in eclipse and out of eclipse, come up with upper limit on tau.
# the graph is then a grid of points of tip and tilt, with upper limits on the mass of the disk in that configuration for a given radius of the disk.

t_mid = 58195. * u.day # time when star has closest approach to planet

r_hill = rhill(bp.M_star, bp.M_c, bp.a_c).to(u.au) # radius of hill sphere

print('Hill sphere radius {}'.format(r_hill))

f_hill = 0.60 # radius of disk in fraction of hill sphere

r_disk = r_hill * f_hill # radius of disk
print('Disk radius {}'.format(r_disk))

angle_i   = 20   # inclination of the disk in degrees
angle_phi = 50   # tilt of the disk in degrees
b = 0.2
impact_b = b * r_hill # impact parameter B of the star

print('Impact parameter distance {}'.format(impact_b))

def ringedge(radius, i_deg, phi_deg, npoints=201):
    'returns evenly spaced points around a ring in the sky coordinates'
    th = np.linspace(0, 2*np.pi, npoints, endpoint=False)
    Xr = radius*np.cos(th)
    Yr = radius*np.sin(th)
    (xs,ys,rs) = ring_to_sky(Xr, Yr, i_deg, phi_deg)
    return (xs,ys)

def xp(t, t_mid, vtrans):
    """Simple function to convert time to x position w.r.t. the planet
        t_mid = u.day, vtrans = u.km/u.s"""
    return ((t-t_mid)*vtrans).to(u.au)


def disk_mass(r_disk, tau, mean_a=0.5*u.micron, mean_rho=2.5*u.g/(u.cm*u.cm*u.cm)):
    'simple mass for a face-on circular optically thin disk'

    # cadged from Mellon's derivation in thesis - p.44, eq. 4.9
    Mdisk = (4*np.pi*mean_a*mean_rho*tau*r_disk*r_disk)/3.
    return Mdisk.to(u.g)


def disk_lc_model(x, xlower, xupper, deltadisk=0.0, foutside=1.0):
    'simple photometric model for a disk'
    indisk = (x>xlower) * (x<xupper)
    f = np.full(x.size, foutside)
    if np.sum(indisk):
        f[indisk] = f[indisk] - deltadisk
    return f
    
n_i =   100
n_phi = 100

i_range = np.linspace(10, 90, n_i)
phi_range = np.linspace(0, 180, n_phi)
phi_min = phi_range.min()
phi_max = phi_range.max()

dphi = 0.5*((phi_max-phi_min)/(phi_range.size-1))

# make a time series dataset
t = np.linspace(t_mid-50*u.day, t_mid+50*u.day, 300)

#determine velocity of planet
vplanet = vcirc(bp.M_star, bp.M_c, bp.a_c)
print('circular velocity at planet is {}'.format(vplanet))

#absorption constant 
taudisk = 0.1


def run_cpd_fit(offset_days):
    #Make mock photometry
    f_sig = 0.01
    f = np.random.normal(1.0, f_sig, t.size)
    ferr = np.full(f.size, f_sig)
    
    # now select the times where the star is within the disk radius
    x_coord = xp(t,t_mid,vplanet)
    (xring, yring, radring) = sky_to_ring(x_coord, impact_b, angle_i, angle_phi)
    
    indisk = (radring < r_disk)
    
    # #Offset photometric data
    t_min = t_mid - (offset_days)*u.day
    t_offset = len(t[(t>t_min) & (t<t_mid)])
    indisk_shift = np.full((t_offset), False, dtype=bool)
    indisk = np.concatenate((indisk[t_offset:], indisk_shift))
    
    # put in an absorption trough
    f[indisk] = f[indisk] - taudisk

    print('estimated thin disk mass: {}'.format(disk_mass(1*u.au, 1e-3)))
    
    # mass estimate from ALMA Perez et al. 2019a.
    print('Mass estimate from lacour2021: {}'.format((bp.M_c * 5e-8).to(u.g)))
    
    gmodel = Model(disk_lc_model, independent_vars=['x','xlower','xupper'])
    print('parameter names: {}'.format(gmodel.param_names))
    print('independent variables: {}'.format(gmodel.independent_vars))
    
    params = gmodel.make_params(deltadisk=0.0, foutside=1.0, xlower=2.5, xupper=3) #, xlower=1.8, xupper=3.4)
    
    deltaf = np.zeros((i_range.size,phi_range.size))
    redchisq = np.zeros_like(deltaf)
    deltaf_err = np.zeros_like(deltaf)
    fout = np.zeros_like(deltaf)
    fout_err = np.zeros_like(deltaf)
    sanity = np.zeros_like(deltaf)
    
    for i, curr_i in enumerate(i_range):
        for j, curr_phi in enumerate(phi_range):
    
            (xring, yring, radring) = sky_to_ring(x_coord, impact_b, curr_i, curr_phi)
            indisk2 = (radring < r_disk)
            nindisk = np.sum(indisk2)
    
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
                sanity[i][j] = curr_i

    # okay, propagate the errors for fout
    tau = deltaf/fout
    tau_err = np.power(np.power((deltaf_err/deltaf),2.)+np.power((fout_err/fout),2),0.5)*tau
    
    f = np.random.normal(1.0, f_sig, t.size)
    f[indisk] = f[indisk] - taudisk
    return deltaf, f

deltaf_all = np.empty([6, n_i, n_phi])
f_all = np.empty([6, t.size])


for i, offset_day in enumerate(offset_days):
    print(f'run CPD model with offset of {offset_day} days from the midpoint at {t_mid}')
    deltaf_all[i,...], f_all[i,...] = run_cpd_fit(offset_day)
    
    
import matplotlib
matplotlib.rcParams.update({'font.size': 16})

def plot_tau(deltaf, f,t, offset_list): 
    t_0 = t_mid  #58000 *u.d
    shrink = 0.8
    fig2, f2_axes = plt.subplots(4,3, figsize=(8.3*1.9, 11.7*0.75), constrained_layout=True, gridspec_kw={'height_ratios':[4,1,4,1]})
    axes = f2_axes.flatten()

    i_min = i_range.min()
    i_max = i_range.max()
    di = 0.5*((i_max-i_min)/(i_range.size-1))
    
    phi_min = phi_range.min()
    phi_max = phi_range.max()
    dphi = 0.5*((phi_max-phi_min)/(phi_range.size-1))
    n=0
    p=0
        
    for ax in axes[:3]:
        ax.set_title(r" $\Delta$t"+f"= {offset_list[n]}d" )
        im1=ax.imshow(deltaf[n],
                 extent=(phi_min-dphi,phi_max+dphi,i_min-di,i_max+di),
                 vmin=0,
                 vmax=taudisk*1.1,
                 origin='lower')
        ax.set_xlabel(r"Tilt $ \phi $ [deg]")
        ax.scatter(50,20, color='red')
        ax.xaxis.set_major_locator(MultipleLocator(30))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_major_locator(MultipleLocator(30))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))

        # For the minor ticks, use no labels; default NullFormatter.
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(10))
        ax.set_xlim(-5,185)
        ax.set_ylim(5,95)
        n += 1
    for ax in axes[6:9]:
        ax.set_title('\n'+r"$\Delta$t"+f"= {offset_list[n]}d" )
        im1=ax.imshow(deltaf[n],
                 extent=(phi_min-dphi,phi_max+dphi,i_min-di,i_max+di),
                 vmin=0,
                 vmax=taudisk*1.1,
                 origin='lower')
        ax.set_xlabel(r"Tilt $ \phi $ [deg]")
        ax.scatter(50,20, color='red')
        ax.xaxis.set_major_locator(MultipleLocator(30))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_major_locator(MultipleLocator(30))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))

        # For the minor ticks, use no labels; default NullFormatter.
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(10))
        ax.set_xlim(-5,185)
        ax.set_ylim(5,95)
        n += 1
    for ax in axes[3:6]:
        ax.axvspan(58193-t_0.value, 58203-t_0.value, alpha=0.3, color= 'orchid')
        ax.axvspan(min(t[f[p]<0.95]).value-t_0.value,max(t[f[p]<0.95]).value-t_0.value, alpha=0.3, color= 'lightskyblue')

        ax.scatter(t-t_0, f[p], s= 10, color='#1f77b4')
        ax.set_xlabel(f"Day [MJD-{int(t_mid.value)}]")
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.xaxis.set_major_locator(MultipleLocator(20))
        p+=1
    for ax in axes[9:]:
        ax.axvspan(58193-t_0.value, 58203-t_0.value, alpha=0.3, color= 'orchid')
        ax.axvspan(min(t[f[p]<0.95]).value-t_0.value,max(t[f[p]<0.95]).value-t_0.value, alpha=0.3, color= 'lightskyblue')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.xaxis.set_major_locator(MultipleLocator(20))
        ax.scatter(t-t_0, f[p], s= 10, color='#1f77b4')
        ax.set_xlabel(f"Day [MJD-{int(t_mid.value)}]")
        p+=1
    cbar = fig2.colorbar(im1, ax=f2_axes, label = r'$\tau\ \sin( \theta)$', shrink=1.0, anchor=(0, 0), pad=0.02) 
    for ax in axes:
        ax.set_aspect('auto')
    axes[0].set_ylabel(r'Inclination $ \theta $ [deg]') 
    axes[6].set_ylabel(r'Inclination $ \theta $ [deg]') 
    axes[3].set_ylabel("Normalized flux ") 
    axes[9].set_ylabel("Normalized flux") 

    plt.draw()
    plt.savefig(paths.figures / 'sensitivity_CPD_mock_data.pdf', bbox_inches='tight')
    
plot_tau(deltaf_all,f_all,t, offset_days)