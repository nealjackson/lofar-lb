import astropy,os,numpy as np
import astropy.io.fits as pyfits

def fits_axzapit (h, k):
    try:
        h.remove(k)
        return False
    except:
        return True

def fits_axzap (infile,outfile):
    hdulist = pyfits.open(infile)
    hdu = hdulist[0]
    while hdu.data.ndim>2:
        shrinkax = np.argwhere(hdu.data.shape==np.min(hdu.data.shape))[0][0]
        hdu.data = np.take (hdu.data,0,axis=shrinkax)
    i=2
    while True:
        i+=1
        isit = fits_axzapit (hdu.header,'CRPIX%d'%i)
        fits_axzapit (hdu.header,'CTYPE%d'%i)
        fits_axzapit (hdu.header,'CUNIT%d'%i)
        fits_axzapit (hdu.header,'CRVAL%d'%i)
        fits_axzapit (hdu.header,'CDELT%d'%i)
        fits_axzapit (hdu.header,'CRPIX%d'%i)
        fits_axzapit (hdu.header,'CROTA%d'%i)
        hdu.header.update()
        if isit:
            break
    os.system ('rm '+outfile)
    hdu.writeto (outfile)

        
