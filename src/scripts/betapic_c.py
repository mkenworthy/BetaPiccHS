import numpy as np
from astropy import constants as c
from astropy import units as u
'contains all constants about beta pic system to be imported and used for all data reduction'
# star
M_star = 1.797 * u.Msun  # Zwintz19
R_star = 1.497 * u.Rsun # Zwintz19
T_star = 8090. * u.K    # Zwintz19
d_star = 19.450 * u.pc # nielsen20
# planet b
# M_b = 12.8 * u.Mjupiter # Nielsen20
# r_b = 1.46 * u.Rjupiter # Chilcote17
# T_b = 1724. * u.K # Chilcote17
# a_b = 10.2 * u.au # Nielsen20
# i_b = 88.82 * u.deg # Nielsen20
# P_b = 22.7 * u.yr # Nielsen20
# d_star = 19.450 * u.pc # nielsen20


# planet c
M_b = 8.5 * u.Mjupiter # Vandal2020
r_b = 1.2 * u.Rjupiter # Nowak2020
T_b = 1250 * u.K # Nowak2020
a_b = 2.68 * u.au # Lacour2021
i_b = 88.95 * u.deg # lacour2021
P_b = 1221 * u.day # lacour2021

"""Calculations for beta pic b"""
# Lagrange (2020)
#
# closest approach 2017 sept 13 = 58009
# 2017 april 11 = 57854
# 2018 feb 16 = 58165

# table containing the Hill sphere transit events for beta pic b
from astropy.table import Table
th = Table()

"""Calculations for beta pic b"""

t = 45/2 #1.5 months 191.8/2
x = 58210
th['th']   = [x-t, x-t/2, x, x+t/2, x+t]
th['therr'] = [18.,   18.,   2.6,   18.,   18.]

th['radius'] = ['100% ingress','50% ingress','mid','50% egress','100% egress']
th.meta['comments'] = ['Beta Pic b Hill sphere transit times',
                       'ingress, 50 percent, closest approach, 50 percent, egress',
                       'times im MJD, errors are on each date']

def addhill(ax, th=th, bottom=-10,height=1000):
    'adds greyed rectangles in Axes with Hill sphere events'
    from matplotlib.patches import Rectangle

    rect100 = Rectangle((th['th'][0], bottom), th['th'][4]-th['th'][0], height,
                            facecolor='0.6',zorder=-6)
    rect50  = Rectangle((th['th'][1], bottom), th['th'][3]-th['th'][1], height,
                            facecolor='0.8',zorder=-5)
    rect0   = Rectangle((th['th'][2]-2, bottom), 4, height,
                            facecolor='0.9',zorder=-4)
    ax.add_patch(rect100)
    ax.add_patch(rect50)
    ax.add_patch(rect0)
    
############################
## TESTING
############################

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    from astropy.time import Time

    t = Time(th['th'], format='mjd', scale='utc')

    print("Hill sphere times")

    for (tn,rn) in zip(t, th['radius']):
        print('{} at {} at a date of {}'.format(rn, tn, tn.iso))

    print(th['th'][4] - th['th'][0])
