#
#  Plots radio sources in a LOTSS field given a file lotss.txt with the following columns:
#
#Source_id,RA,E_RA,E_RA_tot,DEC,E_DEC,E_DEC_tot,Peak_flux,E_Peak_flux,E_Peak_flux_tot,Total_flux,E_Total_flux,E_Total_flux_tot,Resolved,Isl_rms,S_Code,Mosaic_ID

import numpy as np
import astropy,aplpy
import matplotlib
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt, patches as mpatches
import os,sys
plt.rcParams['image.origin'] = 'lower'

def mkfits (rasiz,decsiz,imsiz,pixsiz):
    hdu=pyfits.PrimaryHDU(np.zeros((imsiz,imsiz)))
    hdu.header.update({'CTYPE1':'RA---SIN'})
    hdu.header.update({'CRVAL1':rasiz})
    hdu.header.update({'CRPIX1':imsiz/2.})
    hdu.header.update({'CDELT1':-pixsiz})
    hdu.header.update({'CTYPE2':'DEC--SIN'})
    hdu.header.update({'CRVAL2':decsiz})
    hdu.header.update({'CRPIX2':imsiz/2.})
    hdu.header.update({'CDELT2':pixsiz})
    hdu.header.update({'EQUINOX':2000.0})
    hdu.header.update({'EPOCH':2000.0})
#    hdu.data = np.random.random(imsiz*imsiz).reshape(imsiz,imsiz)
    os.system('rm temp.fits')
    hdu.writeto('temp.fits')


LABELS,RA,DEC,PFLUX,TFLUX,CODE = 0,1,4,7,10,15
a = np.loadtxt('lotss.txt',skiprows=1,delimiter=',',dtype='S')
ra = np.asarray(a[:,RA],dtype='float')
dec = np.asarray(a[:,DEC],dtype='float')
pflux = np.asarray(a[:,PFLUX],dtype='float')
tflux = np.asarray(a[:,TFLUX],dtype='float')
labels = np.asarray(a[:,LABELS])
code = np.asarray(a[:,CODE])
ra=ra[code=='S']
dec=dec[code=='S']
pflux=pflux[code=='S']
tflux=tflux[code=='S']
labels=labels[code=='S']
dec_m, dec_r = dec.mean(), dec.max()-dec.min()
ra_m, ra_r = ra.mean(), (ra.max()-ra.min())*np.cos(np.deg2rad(dec_m))
isize = 1.1*max(dec_r,ra_r)
im = mkfits (ra_m,dec_m, 256, isize/256.0)
gc = aplpy.FITSFigure('temp.fits')
gc.add_grid()
gc.grid.set_color('black')
gc.show_circles(ra,dec,0.03*np.log10(pflux),color='red')
for i in range (min(50,len(ra))):
    astr = labels[i][4:8]+labels[i][12:17]
    xy = [ra[i]+1e-5,dec[i]+1e-5]
    gc.add_label (ra[i]-0.35,dec[i]+0.05,astr,size=7)

gc.save('lotss_plot.png')
