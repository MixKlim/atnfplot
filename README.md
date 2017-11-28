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
