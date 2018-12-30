#
#   aplpy_makeplot.py
#
#   Uses sextractor to find components and then makes a plot which includes them
#
#   arguments:
#   1) name of the fits file
#   2) docut - controls the sigma-level for finding the components. Default is 2.
#              >2 finds fewer components (higher cut), <2 more. <0 probably not useful.
#              If <-1.0, is interpreted as a request for a radius of abs(docut) arcsec.
#   3) outpng - the output png, default is infile with fits replaced by png
#   4) nolabel - do not add label with rms etc
#   5) crms - lowest contour level as a multiple of rms, default 3
#   6) noshift - do not shift centre of plot from centre of image
import numpy as np,matplotlib,scipy
from matplotlib import pyplot as plt; from scipy import stats
import astropy,aplpy,sys,os,argparse
import astropy.io.fits as pyfits
import fits_axzap
sextractor = '/home/jackson/sextractor/usr/bin/sex'
f=open('default.param','w')
f.write('NUMBER\nX_IMAGE\nY_IMAGE\nFLUX_ISO\nFLUX_RADIUS\nFLUX_MAX\nISOAREA_IMAGE\n')
f.close()
os.system(sextractor+' -d >default.sex')
f=open('default.conv','w')
f.write('CONV NORM\n# 3x3 ``all-ground\'\' convolution mask with FWHM = 2 pixels\n')
f.write('1 2 1\n2 4 2\n1 2 1\n')
f.close()

parser = argparse.ArgumentParser()
parser.add_argument ('infits')#
parser.add_argument ('-c', '--docut',type=float,help='higher values include only more significant components; if <-1.0, radius is abs(docut) arcsec',default=2.0)
parser.add_argument ('-o', '--outpng',help='name of output png file',default='')
parser.add_argument ('-r', '--crms',help='multiple of rms for lowest contour',type=float,default=3.0)
parser.add_argument ('-l', '--nolabel',dest='nolabel',help='do not write useful information on plot',action='store_true')
parser.add_argument ('-s', '--noshift',dest='noshift',help='do not shift centre',action='store_true')
parser.add_argument ('-m', '--margin',type=float,help='margin allowance',default=1.7)
args = parser.parse_args()

def make_cut_plot (infits, docut=2.0, outpng='', nolabel=False,  crms=3.0, noshift=False, margin=1.7):
    print 'Plotting ',infits
    fits_axzap.fits_axzap (infits,'temp.fits')
    a = pyfits.getdata('temp.fits')
    h = pyfits.getheader('temp.fits')
    nx,ny = h['NAXIS1'],h['NAXIS2']
    field_radius = h['CDELT2']*ny/2.0
    nx,ny = a.shape[1],a.shape[0]
    trms,tmax = np.array([]), np.array([])
    for i in range(10):
        x1,x2 = int(0.1*i*nx),int(0.1*(i+1)*nx-1)
        for j in range(10):
            y1,y2 = int(0.1*j*ny),int(0.1*(j+1)*ny-1)
            trms = np.append(trms,np.std(a[y1:y2,x1:x2]))
            tmax = np.append(tmax,np.std(a[y1:y2,x1:x2]))
    rms = np.nanmedian(trms)
    vmin,vmax = np.nanmin(a),np.nanmax(a)
    pyfits.writeto('temp.fits',a,h,overwrite=True)
    os.system(sextractor+' temp.fits')
    s = np.loadtxt('test.cat')
    s = np.array([s]) if s.ndim==1 else s
    reqsig = stats.norm.ppf(1-0.5/float(nx*ny)) + (2.0 if docut<-1.0 else docut)
    s = s[s[:,5]>reqsig*rms]
    for i in s:
        print 'Component at (%.1f,%.1f): %f' % (i[1],i[2],i[5])
    gc = aplpy.FITSFigure('temp.fits')
    pixra,pixdec = np.mean(s[:,1]),np.mean(s[:,2])
    dec = h['CRVAL2']+h['CDELT2']*(pixdec-h['CRPIX2'])
    cosdec = np.cos(np.deg2rad(dec))
    ra = h['CRVAL1']+h['CDELT1']*(pixra-h['CRPIX1'])/cosdec
    deccen = h['CRVAL2']+h['CDELT2']*(0.5*ny-h['CRPIX2'])
    cosdeccen = np.cos(np.deg2rad(deccen))
    racen = h['CRVAL1']+h['CDELT1']*(0.5*nx-h['CRPIX1'])/cosdeccen
    pix_range = max(s[:,1].max()-s[:,1].min(),s[:,2].max()-s[:,2].min())
    deg_range = max(margin*h['CDELT2']*pix_range,0.1*ny*h['CDELT2'])
    if docut<-1.0:
        deg_range = -2.0*docut/3600.
    print 'Range is %.1f arcsec, %.1f pix\n'%(1800.0*deg_range,pix_range)
    if noshift:
        gc.recenter(racen,deccen,0.5*deg_range)
    else:
        gc.recenter(ra,dec,0.5*deg_range)
    gc.set_tick_color('black')
    gc.show_colorscale(cmap=matplotlib.cm.gray_r,vmin=vmin,vmax=vmax)
    levels, tlevels = [vmax], vmax
    while tlevels > crms*rms:
        tlevels /= np.sqrt(2)
        levels.append(tlevels)
    gc.show_contour(levels=np.sort(levels))
    if not nolabel:
        bstr=''
        for i in range(len(h)):
            try:
                if 'BMAJ' in h[i] and 'AIPS' in h[i]:
                    bmaj = 3600.*float(h[i].split('BMAJ=')[1].split()[0])
                    bmin = 3600.*float(h[i].split('BMIN=')[1].split()[0])
                    bstr = 'beam %.1fx%.1f'%(bmaj,bmin)
            except:
                pass
        if rms>0.1:
            gc.add_label(0.5,0.05,'Peak %.1f, rms %.1f Jy %s'%\
                    (vmax,rms,bstr),relative=True,size=14)
        elif rms>1.e-4:
            gc.add_label(0.5,0.05,'Peak %.1f, rms %.1f mJy %s'%\
                    (vmax*1000.,rms*1000.,bstr),relative=True,size=14)
        else:
            gc.add_label(0.5,0.05,'Peak %.1f, rms %.1f uJy %s'%\
                    (vmax*1.e6,rms*1.e6,bstr),relative=True,size=14)
    if outpng=='':
        outpng=infits.replace('fits','png')
    gc.save(outpng)
    os.system('rm default.conv;rm default.sex;rm default.param;rm temp.fits')
    

make_cut_plot (args.infits, args.docut, args.outpng, args.nolabel, args.crms, args.noshift, args.margin)
