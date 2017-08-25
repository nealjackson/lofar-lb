# spflg  Neal Jackson 25/8/17
import sys,os,idx_tels,numpy as np, matplotlib, time
import pyrap; from pyrap import tables as pt
import matplotlib;from matplotlib import pyplot as plt

def dget_t (vis, tel1, tel2, datacolumn, spwsel=np.nan):
    os.system('taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' % (vis, tel1, tel2, 'cl_temp.ms'))
    t = pt.table('cl_temp.ms')
    dc = t.select(datacolumn)
    d = np.asarray([tuple(each.values()) for each in dc])[:,0,:,:]
    d = np.swapaxes(d,0,2)     # makes the order pol - chan - time as in casa
    spw = np.ravel(np.asarray([tuple(each.values()) for each in \
              t.select('DATA_DESC_ID')]))
    if ~np.isnan(spwsel) and spwsel in np.unique(spw):    # request just one spw
        return d[:,:,spw==spwsel]
    if spw.sum():
        for i in np.unique(spw):
            new = np.take(d,np.argwhere(spw==i),axis=2)[:,:,:,0]
            try:
                d_out = np.concatenate((d_out,new),axis=1)
            except:
                d_out = np.copy(new)
        d = d_out
    return d

def spflg (vis,tel=['ST001','DE601'],lastv=-1,datacolumn='DATA',pol=0,doplot=False):
    # Find which number is which antenna
    itels = np.sort(idx_tels.get_idx_tels (vis, tel))
    if itels == []:
        print 'Failed to find telescopes, returning'
        return -1,-1
    d1 = dget_t (vis,itels[0],itels[1],datacolumn=datacolumn)
    d1=d1[pol,:,:lastv]
    a1,p1 = np.abs(d1),np.arctan2(d1.imag,d1.real)
    a2 = np.sort(np.ravel(a1))
    amax = a2[99*len(a2)/100]
    if doplot:
        plt.subplot(121)
        plt.imshow(a1.T,vmax=amax,origin='lower',aspect='auto',\
                   interpolation='none',cmap=matplotlib.cm.hot)
        plt.title(vis+' '+tel[0]+'-'+tel[1])
        plt.xlabel('Amplitude')
        plt.subplot(122)
        plt.imshow(p1.T,origin='lower',aspect='auto',interpolation='none',\
                   cmap=matplotlib.cm.hot)
        plt.xlabel('Phase')
        plt.show()
    else:
        return a1.T, p1.T
