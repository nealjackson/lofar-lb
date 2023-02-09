#
#  Plots a field, showing WENSS sources and LBCS calibrators
#  LBCS calibrators are plotted larger and redder the more antennas have
#  coherent phase solutions (i.e. large red ones can be used for all baselines)
#  Wenss sources are plotted according to flux
#  If LBCS and WENSS catalogues are not present they will be downloaded.
#  Invoke: python lbcs_plot.py [ra-centre] [dec-centre]
import numpy as np
import astropy,aplpy
import matplotlib
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt, patches as mpatches
from correlate import *
import os,sys
plt.rcParams['image.origin'] = 'lower'
BPIX = 600

def mkfits (rasiz,decsiz,imsiz,pixsiz):
    hdu=pyfits.PrimaryHDU(np.zeros((int(imsiz),int(imsiz))))
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


try:
    lbcs = np.loadtxt('lbcs_stats.sum',dtype='S')
except:
    os.system('wget http://www.jb.man.ac.uk/~njj/lbcs_stats.sum')
    lbcs = np.loadtxt('lbcs_stats.sum',dtype='S')

try:
    w = np.loadtxt('wenss2000.simple',dtype='S')
except:
    os.system('wget http://www.jb.man.ac.uk/~njj/wenss2000.simple')
    w = np.loadtxt('wenss2000.simple',dtype='S')

w = w[np.asarray(w[:,9],'int')>99]
wra = 15.*np.asarray(w[:,1],dtype='f')+0.25*np.asarray(w[:,2],dtype='f')+\
        np.asarray(w[:,3],dtype='f')/240.
wdec = np.asarray(w[:,4],dtype='f')+np.asarray(w[:,5],dtype='f')/60.+\
        np.asarray(w[:,6],dtype='f')/3600.
wf = np.asarray(w[:,9],dtype='f')
wc = np.dstack((wra,wdec))[0]
coords = np.asarray(lbcs[:,-2:],dtype='float')
ra,dec = float(sys.argv[1]),float(sys.argv[2])
a = correlate (np.array([[ra,dec]]),0,1,coords,0,1,3.0)
wa = correlate (np.array([[ra,dec]]),0,1,wc,0,1,3.0)
l = np.take(lbcs,np.asarray(a[:,1],dtype='int'),axis=0)
wl = np.take(w,np.asarray(wa[:,1],dtype='int'),axis=0)
wra = np.take(wra,np.asarray(wa[:,1],dtype='int'),axis=0)
wdec = np.take(wdec,np.asarray(wa[:,1],dtype='int'),axis=0)
wf = np.take(wf,np.asarray(wa[:,1],dtype='int'),axis=0)
mkfits (ra, dec, BPIX, 1.0/100)
lra,ldec = np.asarray(l[:,7],dtype='f'),np.asarray(l[:,8],dtype='f')
lnames = l[:,0]
coh = np.array([])
for i in l[:,5]:
    coh = np.append(coh, 1.0*i.count('P')/\
                    (i.count('P')+i.count('S')+i.count('X')))

gc = aplpy.FITSFigure('temp.fits')
gc.add_grid()
gc.grid.set_color('black')
for i in range(1,4):
    gc.show_circles (ra,dec,i,color='red')

for i in range(len(wf)):
    gc.show_circles(wra[i],wdec[i],0.05,color='green',facecolor='green',alpha=min(0.001*wf[i],1.0))

for i in range(len(l)):
    if coh[i]>0.0:
        gc.show_circles(lra[i],ldec[i],0.1*coh[i],color='red',facecolor='red',alpha=coh[i])
        gc.add_label(lra[i]-0.5,ldec[i]+0.15,lnames[i])


        
gc.save('lbcs_plot.png')
os.system('display lbcs_plot.png')
print 'Please make a note of the L-numbers you want to use, and then close the plot'
print 'Now enter the L-numbers one by one (e.g. \'L444444\'); X to finish'
f=open('lbcs_cat.csv','w')
f.write('raj2000,decj2000,ObsID\n')
while True:
    s = raw_input('enter source> ')
    if not len(s):
        continue
    if s[0]=='X' or s[0]=='x':
        break
    try:
        idx = np.argwhere(l[:,0]==s)[0][0]
        print '%f,%f,%s\n'%(lra[idx],ldec[idx],s)
        f.write('%f,%f,%s\n'%(lra[idx],ldec[idx],s))
    except:
        print 'could not find ',s
        print l[:,0]

f.close()
print 'Written lbcs_cat.csv'
