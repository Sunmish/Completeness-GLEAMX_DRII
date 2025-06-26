#!/usr/bin/env python

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import math as m
from optparse import OptionParser




# Read input parameters
usage="Usage: %prog [options] <srcname> <map> <map_psf> <output_file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--z',type="float", dest="z", default=-26.7, help="Dec at zenith, in deg [default=%default]")
parser.add_option('--flux',type="float", dest="flux", default=1.0, help="Flux density of simulated sources, in Jy [default=%default]")
(options, args) = parser.parse_args()
srcname=args[0]; map=args[1]; map_psf=args[2]; output_file=args[3]
z=options.z
flux=options.flux

# Define output file
f=open(output_file,'w')

print('# RA Dec S bmaj bmin bpa R', file=f)

# -------------------------------------------

# Read file with list of source positions in deg
ra_list=[]
dec_list=[]
for line in open(srcname):
  columns = line.split()
  if columns[0] != '#':
    ra_list.append(float(columns[0]))
    dec_list.append(float(columns[1]))

# Read bmaj and bmin from fits header of surface brightness map
img = fits.open(map)
bmaj=img[0].header['BMAJ']
bmin=img[0].header['BMIN']

# Convert bmaj and bmin from deg to arcsec, and print to screen
bmaj=bmaj*3600.0
bmin=bmin*3600.0
print('bmaj / arcsec =',bmaj)
print('bmin / arcsec =',bmin)

# Read PSF map
img_psf = fits.open(map_psf)
mappix = (np.squeeze(img_psf[0].data))
img_psf.close()
header_psf = fits.getheader(map_psf)
psf_data = fits.getdata(map_psf)
wcs = WCS(header_psf).celestial

for (ra,dec) in np.nditer((ra_list,dec_list)):  
    
    y, x = wcs.all_world2pix(ra, dec, 0)

    y = int(y)
    x = int(x)

    max_x = psf_data[0, ...].shape[0]
    max_y = psf_data[0, ...].shape[1]

    if x >= max_x:
        x_i = max_x-1
    else:
        x_i = x
    if y >= max_y:
        y_i = max_y-1
    else:
        y_i = y

    a = psf_data[0, x_i, y_i]*3600.
    b = psf_data[1, x_i, y_i]*3600.
    c = psf_data[2, x_i, y_i]
    r = 1.
    print(ra,dec,flux,a,b,c,r, file=f)

    
