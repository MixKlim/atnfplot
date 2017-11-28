__author__ = 'mixklim'

"""
The script constructs a beam mosaic for LOFAR observations based on the observational settings 

Input:
01. Observation ID
02. A FITS image file of the observable object in RA / DEC coordinates
03. Coordinates of the observed stellar object [columns: RA in deg, DEC in deg, beam number]
04. A parset file with all necessary information about a given LOFAR observation (central beam coordinates, number of stations, see below)

Output:
A .pdf file of a RA / DEC observational plane tiled by LOFAR beams

The script uses beam geometry calculations from plot_LOFAR_TA_multibeam.py by Joeri van Leeuwen and can be easily adapted to any other telescope

The script was used in the following papers:

Mikhailov, K., van Leeuwen, J., 2016, A&A, 593, 21
The LOFAR search for radio pulsars and fast transients in M 33, M 81, and M 82

van Leeuwen J., Mikhailov K., Keane E., Kondratiev V., Michilli D., Hessels J., LOFAR Pulsar PWG, 2017, A&A, in prep.
Andromeda and the Crab – A LOFAR radio search for bright single pulses and periodic signals from M 31

"""

import aplpy
import astropy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import matplotlib.ticker as ticker
from matplotlib.patches import Ellipse, Rectangle
import numpy as np
import sys, re, os
import ephem
from math import cos, atan, pi
from operator import indexOf
from optparse import OptionParser
from astropy import log
log.setLevel('WARNING')
log.setLevel('INFO')

# Set your input parameters here
obsid = 'L??????' # LOFAR OBS ID
fitsfile = '??????.fits' # FITS file of the observable object in RA / DEC coordinates
datafile = '??????.txt' # a text file with proper beam coordinates
imagefile = 'L246341.pdf' # output PDF image
parset = 'L246341.parset' # LOFAR parset file with required information about observational settings

# load FITS image using APLpy (https://aplpy.github.io/)
gc = aplpy.FITSFigure(fitsfile)
gc.show_grayscale()
gc.add_grid()

##### Add a beam shape to the current figure #####

RAD2DEG = 1. / 0.0174532925
sap = 0  # number of SAPs (LOFAR subarray pointings)
nisbeams = 0 # also check for IS beams

# Get observational settings from the LOFAR parset
file = open(parset, "r")
for line in file.readlines():
    # check "coherent" attribute of the beam
    if line.find("Observation.Beam[%s].Tied" % (sap)) >= 0 and line.find("coherent") >= 0:
        beamtype = line.split()[-1].lower()[:1]
        if beamtype == 'f': nisbeams += 1
    # get the central coordinate
    if line.find("Observation.Beam[%s].angle1" % (sap)) >= 0:
        rarad = np.float64(line.split("=")[-1])
    if line.find("Observation.Beam[%s].angle2" % (sap)) >= 0:
        decrad = np.float64(line.split("=")[-1])
    # get some numbers for the beam estimate
    if line.find("OLAP.storageStationNames=[") >= 0:
        stationlist = line.split("=")[-1]
    if line.find("Observation.startTime") >= 0:
        # pyephem needs yyyy/mm/dd format
        # correct time format is e.g. starttime = '2015/10/15 04:00:00'
        starttime = line.split("=")[-1].replace('-', '/')
    if line.find("Observation.Beam[%s].subbandList" % (sap)) >= 0:
        subbandlist = line.split("=")[-1]
file.close()

# convert to HHMMSS
src = ephem.FixedBody()
src._ra = float(rarad)
src._dec = float(decrad)

# Estimate the beam shape, only for HBA data
# set up the observatory
src._epoch = ephem.J2000
telescope = ephem.Observer()
telescope.lat = ephem.degrees('52.91511')
telescope.long = ephem.degrees('6.869883')
telescope.elevation = 0
telescope.date = ephem.Date(starttime)
src.compute(telescope)

# calculate the approximate beam size
LOW_FREQ = ????? 		# SET YOUR BOTTOM OBS FREQUENCY HERE!
BW = ????? 			# SET YOUR BANDWIDTH HERE!

# CHAN_BW = ?????		# SET YOUR CHANNEL BANDWIDTH HERE!
#BOTTOMSUBBAND = float(re.findall(r"\d+", subbandlist)[0])
#l = 300.0 / (LOW_FREQ + chan_bw * BOTTOMSUBBAND) # bottom wavelength in m

l = 300.0 / (LOW_FREQ + BW / 2) # high wavelength in m 
diameter = 300 # superterp [https://www.astron.nl/radio-observatory/astronomers/users/technical-information/lofar-array-configuration/lofar-array-conf]
beamsize = 1.22 * l / diameter  # FWHM in rad

# calculate the beam angle in ra-dec plot
raoff = telescope.radec_of(src.az, src.alt + beamsize)[0] - telescope.radec_of(src.az, src.alt)[0]
decoff = telescope.radec_of(src.az, src.alt + beamsize)[1] - telescope.radec_of(src.az, src.alt)[1]

beamsize = beamsize * RAD2DEG

# plot the beam shape in the background
angle = atan(raoff / decoff)

gc.add_beam(major=beamsize, minor=beamsize/cos(pi/2-src.alt), angle=-angle * RAD2DEG, corner='bottom left', frame=False,
            borderpad=0.5, pad=0.5)

# Set the transparency level of the beam
gc.beam.set_alpha(0.3)

# Set the color of the whole beam, or the edge and face color individually
gc.beam.set_edgecolor('white')
gc.beam.set_facecolor('green')

# add LOFAR beams
data = np.loadtxt(datafile)
ra, dec, ind = data[:, 0], data[:, 1], data[:, 2]
gc.show_ellipses(ra, dec, height=beamsize, width=beamsize/cos(pi/2-src.alt), angle=-angle * RAD2DEG, layer='scatter_set_1', edgecolor='white', facecolor='none', alpha=0.6, zorder=None)

# add beam numbers
for i in range(0,len(ind)):
	gc.add_label(ra[i], dec[i], int(ind[i]), relative=False, color='orange', family='serif', style='normal', variant='normal', stretch='normal', weight='medium', size='x-large', fontproperties=None, horizontalalignment='center', verticalalignment='center', layer=None)

# font properties
gc.axis_labels.set_font(size='large', weight='bold', stretch='normal', family='sans-serif', style='normal', variant='normal')
             
gc.save(imagefile)
