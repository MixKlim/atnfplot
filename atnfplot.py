__author__ = 'mixklim'

"""
The script constructs a P-Pdot diagram of known pulsars, RRATs, and magnetars

Input:
01. An updated link to the ATNF pulsar catalogue
02. An updated link to the RRATalog
03. An updated link to the McGill Online Magnetar Catalogue

Output:
A pdf plot of the diagram

The script was used in the following papers:

Mikhailov, K., van Leeuwen, J., Jonker, P. G., 2017, ApJ, 840, 9
A Search for Millisecond­pulsar Radio Emission from the Faint Quiescent Soft X­Ray Transient 1H 1905+000

"""

import urllib2
import numpy as np
import pylab as pl

####################################################
#		   Normal Pulsars 		   #
####################################################

url1 = 'http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?version=1.57&Name=Name&P0=P0&P1=P1&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names=&ephemeris=long&submit_ephemeris=Get+Ephemeris&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+last+digit+error&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query'

data1 = urllib2.urlopen(url1)

np0 = []
np1 = []

# skip first 33 lines with meta-data
for i in xrange(33):
	data1.next()

# read all the relevant normal pulsars from the ATNF catalogue
for line in data1:
    if "-----" not in line:
    	if not line.isspace():
    		words = line.split()
    		if words[5] is not '*' and words[9] is not '*':
    			np0.append(words[5])
    			np1.append(words[9])
    else: break

#################################################
#			RRATs 			#
#################################################

url2 = 'http://astro.phys.wvu.edu/rratalog/rratalog.txt'

data2 = urllib2.urlopen(url2)

rr0 = []
rr1 = []

# skip first line with meta-data
for i in xrange(1):
	data2.next()

# read all the relevant RRATs from the rratalog
for line in data2:
	words = line.split()
    	if '-' not in words[1] and '-' not in words[2]:
    		rr0.append(words[1])
    		rr1.append(words[2]+'E-15')

#################################################
#		     Magnetars 			#
#################################################

url3 = 'http://www.physics.mcgill.ca/~pulsar/magnetar/main.html'

from bs4 import BeautifulSoup

data3 = urllib2.urlopen(url3)
soup = BeautifulSoup(data3)

lines = []
mg0 = []
mg1 = []

for tag in soup.find_all('tr'):
        lines.append(tag.get_text())

# read all the relevant magnetars from the McGill Online Magnetar Catalogue
for line in lines:
    words = line.split()
    if words: # skip empty lines
        if 'Name' not in words and '(J2000)' not in words: # skip header lines
            if '...' not in words[2] and '...' not in words[4]:
                mg0.append(str(words[2]).translate(None, '[<~!()*]'))
                mg1.append(str(words[4]).translate(None, '[<~!()*]')+'E-11')

################################################
#                make pulsar plots             #
################################################

plot1, = pl.loglog(np0, np1, 'ro', markersize=2, label='Radio Pulsar')
plot2, = pl.loglog(rr0, rr1, 'y^', markersize=8, label='RRAT')
plot3, = pl.loglog(mg0, mg1, 'g*', markersize=10, label="Magnetars")

# log-scale grid
x = np.logspace(-3.0, np.log10(20), num=50, endpoint=True, base=10.0)

################################################
#	      add magnetic field B 	       #
################################################

def field(B, x):
    return pow(B, 2) / (x * pow(3.2e19, 2))
# make magnetic field plots
B = 10**8
while B < 10**16:
    plot, = pl.loglog(x, field(B, x), 'k--')
    B *= 10

################################################
#	     add characteristic age 	       #
################################################

def age(t, x):
    return x / (2 * t)
# make age plots
T0 = 3600 * 24 * 365.25 * 1000 * 100
T = T0
while T < T0 * 1e5:
    plot, = pl.loglog(x, age(T, x), 'k:')
    T *= 100

#####################################################################################
#		add pulsar death line [eq .(9) from Chen & Ruderman 1993] 	    #
#####################################################################################

C = 3.2e19
plot, = pl.loglog(x, 10**((78 - 7 * np.log10(C) + 9.5 * np.log10(x)) / 3.5), 'k-')

#################################################
#                    plot design                #
#################################################

pl.rc('text', usetex=True)
pl.rc('font', family='serif')

# give plot a title
pl.title('$P-\dot{P}$ diagram', fontsize=16)

# make axis labels
pl.xlabel(r'Period (s)', fontsize=16)
pl.ylabel(r'Period Derivative (s s$^{-1}$)', fontsize=16)

# set axis limits
pl.xlim(10**(-3), 20)
pl.ylim(10**(-21), 10**(-9))

# set axis ticks size
pl.xticks(fontsize=15)
pl.yticks(fontsize=15)

#################################################
#                    text labels                #
#################################################

# magnetic fields
pl.text(pow(10,0.05), pow(10,-9.5),r'$\bf{10^{15}}~\mbox{G}$',rotation=-15)
pl.text(pow(10,-1), pow(10,-10.47),r'$\bf{10^{14}}~\mbox{G}$',rotation=-15)
pl.text(pow(10,-2.85), pow(10,-10.65),r'$\bf{10^{13}}~\mbox{G}$',rotation=-15)
pl.text(pow(10,-2.85), pow(10,-12.65),r'$\bf{10^{12}}~\mbox{G}$',rotation=-15)
pl.text(pow(10,-2.85), pow(10,-14.65),r'$\bf{10^{11}}~\mbox{G}$',rotation=-15)
pl.text(pow(10,-2.85), pow(10,-16.65),r'$\bf{10^{10}}~\mbox{G}$',rotation=-15)
pl.text(pow(10,-2.85), pow(10,-18.65),r'$\bf{10^{9}}~\mbox{G}$',rotation=-15)
pl.text(pow(10,-2.85), pow(10,-20.6),r'$\bf{10^{8}}~\mbox{G}$',rotation=-15)

# characteristic ages
pl.text(pow(10,-1.9), pow(10,-14.22),r'$\bf{100}~\mbox{kyr}$',rotation=17)
pl.text(pow(10,-1.9), pow(10,-16.23),r'$\bf{10}~\mbox{Myr}$',rotation=17)
pl.text(pow(10,-1.9), pow(10,-18.3),r'$\bf{1}~\mbox{Gyr}$',rotation=17)

# death line
pl.text(pow(10,-1.05), pow(10,-18.8),r'\textit{Death~line}',rotation=36)

#################################################
#                    plot legend               #
#################################################

# make legend
legend = pl.legend(loc='lower right', shadow=True, numpoints=1)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

# save plot (pdf format)
pl.savefig('PPdot', papertype = 'a4', orientation = 'portrait', format = 'pdf')
