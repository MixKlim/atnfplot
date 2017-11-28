__author__ = 'mixklim'

"""
The script constructs an aitoff map of all known FRBs and RRATs to date. The DM of each source is represented by the size of the point, the associated flux density is represented by the color scheme.

Input:
01. Type of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
02. An updated .csv table of all known FRB (can be found on the FRB Online Catalogue: http://www.frbcat.org/)

Output:
A pdf image of an FRB / RRAT map

The script was used in the following papers:

van Leeuwen J., Mikhailov K., Keane E., Kondratiev V., Michilli D., Hessels J., LOFAR Pulsar PWG, 2017, A&A, in prep.
Andromeda and the Crab – A LOFAR radio search for bright single pulses and periodic signals from M 31
 
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
import ephem
import matplotlib as mpl
import urllib2

# first, choose a type of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
projection = 'aitoff' 
fig = plt.figure(figsize=(11.69,8.27)) # landscape scale
ax = fig.add_subplot(111, projection=projection, axisbg='White')

# create a new cmap with limited color range
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

# set the color maps for different types of sources
frb_cmap = truncate_colormap(plt.cm.get_cmap('YlOrBr'), 0.2, 1.0)
rrat_cmap = truncate_colormap(plt.cm.get_cmap('Reds'), 0.2, 1.0)

#################################################
#                    plot design                #
#################################################

# activate latex text rendering
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=2)
plt.rc('font', family='serif', weight='bold')
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}']

tick_labels = np.array(['10h','8h','6h','4h','2h','0h','22h','20h','18h','16h','14h'])
#tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
#tick_labels = np.remainder(tick_labels + 360 + org, 360)
ax.set_xticklabels(tick_labels)  # we add the scale on the x axis
ax.set_xlabel(r'\bf{Ra}', fontsize=12)
ax.set_ylabel(r'\bf{Dec}', fontsize=12)
ax.grid(True)

#################################################################################
#   scatter plot of points on the celestial sphere using Mollweide projection   #
#################################################################################

def plot_mwd(ra,dec,org=0, color=None, vmin=0, vmax=1, colormap='Greys', shape=10, edgecolor=None, flag=False, fraction=None, pad=None, text=None):
    ''' ra, dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    '''
    x = np.remainder(ra+360-org,360) # shift RA values
    ind = x > 180
    x[ind] -= 360    	# scale conversion to [-180, 180]
    x = -x    		# reverse the scale: East to the left
    # convert degrees to radians
    sc = ax.scatter(np.radians(x), np.radians(dec), c=color, vmin=vmin, vmax=vmax, cmap=colormap, s=shape, edgecolors=edgecolor)  

    # add colormap
    if flag is True:
        c = plt.colorbar(sc, orientation='horizontal', fraction=fraction, pad=pad)
        c.set_label(text, size=15)

#################################################################################
"""
The plane of our galaxy, the Milky Way, in a plot of the celestial sphere in equatorial coordinates RA, Dec
The galactic plane has latitude coordinate b = 0 in the galactic coordinate system.
We will generate an array of galactic coordinates with 0 latitude
and we'll convert it to equatorial coordinates using the Ephem library.
"""
#################################################################################

lon_array = np.arange(0,360)
lat = 0.
eq_array = np.zeros((360,2))
for lon in lon_array:
    ga = ephem.Galactic(np.radians(lon), np.radians(lat))
    eq = ephem.Equatorial(ga)
    eq_array[lon] = np.degrees(eq.get())
Ra = eq_array[:,0]
Dec = eq_array[:,1]

plot_mwd(Ra, -Dec, org=180, color='b')

#################################################################################
#                                   Add FRBs                                    #
#################################################################################

CSV = "frbcat_2017-11-22.csv" # check the latest version on http://www.frbcat.org/
data = ascii.read(CSV, header_start=0, data_start=1)

'''
Plot FRB data. RA, DEC are arrays of the same length. RA takes values in [0,360), Dec in [-90,90],
which represent angles in degrees.
'''

FRB = data["Name"]
frb_flux = data["Flux"]
frb_dm = data["DM"]

ra = coord.Angle(data["RAJ"], unit=u.hour)
RA = ra.degree
dec = coord.Angle(data["DECJ"], unit=u.degree)
DEC = dec.degree

plot_mwd(RA, DEC, color=frb_flux, colormap=frb_cmap, vmin=min(frb_flux), vmax=max(frb_flux), shape=frb_dm, edgecolor='black', flag=True, fraction=0.1, pad = 0.06, text="FRBs flux densities (Jy)")

#################################################################################
#                                   add RRATs                                   #
#################################################################################

url = 'http://astro.phys.wvu.edu/rratalog/rratalog.txt'
data = urllib2.urlopen(url)

'''
Plot RRAT data. RA, DEC are arrays of the same length. RA takes values in [0,360), Dec in [-90,90],
which represent angles in degrees.
'''

rrat_name = []
rrat_dm = []
rrat_ra = []
rrat_dec = []
rrat_flux = []

# skip first line with meta-data
for i in xrange(1):
    data.next()

# read all the relevant RRATs from the rratalog
for line in data:
    words = line.split()
    rrat_name.append(words[0])
    rrat_dm.append(float(words[3]))
    rrat_ra.append(words[4])
    rrat_dec.append(words[5])
    if '-' not in words[12]:
        rrat_flux.append(float(words[12]) * 1e-3) # fluxes in Jy
    else:
        rrat_flux.append(0)

# fix broken DEC elements, e.g. 06:
for dec in rrat_dec:
    if dec[len(dec)-1] is ':':
        rrat_dec[rrat_dec.index(dec)] = dec[:-1]

ra = coord.Angle(rrat_ra, unit=u.hour)
RA = ra.degree
dec = coord.Angle(rrat_dec, unit=u.degree)
DEC = dec.degree

plot_mwd(RA, DEC, color=rrat_flux, colormap=rrat_cmap, vmin=min(rrat_flux), vmax=max(rrat_flux), shape=rrat_dm, edgecolor='black', flag=True, fraction=0.1, pad = 0.04, text="RRATs flux densities (Jy)")

#################################################################################
#                                   the legend                                  #
#################################################################################

# we combine both FRB and RRATS DMs
frb_dm = [dm for dm in frb_dm]

total_dm = rrat_dm + frb_dm
# here we step over the sorted data into 15 strides and select the
# last 5 steps as a representative sample, this only works if your
# data is fairly uniformly distributed
legend_sizes = np.sort(total_dm)[::len(total_dm) // 15][-5:]
# get the indices for each of the legend sizes
indices = [np.where(total_dm == v)[0][0] for v in legend_sizes]
# plot each point again, and its value as a label
for i in indices:
    ax.scatter(100, 105, s=total_dm[i], c='white', label='{:.2f} pc/cc'.format(round(total_dm[i], -2)))

# make legend
legend = ax.legend(loc=(0.0, 1.1), ncol=len(indices), shadow=True, fancybox=True, prop={'size': 25}, scatterpoints=1, handletextpad=0.6, handlelength=0.2)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.8')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('medium')

# The legend line width
for label in legend.get_lines():
    label.set_linewidth(1.5)  

#################################################################################
#                                   final plot                                  #
#################################################################################
plt.savefig("frbmap.pdf", papertype = 'a4', orientation = 'portrait', format = 'pdf')
