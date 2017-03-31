import sys,os,idx_tels,numpy as np, matplotlib, time
from matplotlib import pyplot as plt
try:
    import pyrap; from pyrap import tables as pt
    have_tables = True
except:
    have_tables = False
    import casac, string, inspect
    def define_tb():
        a=inspect.stack()
        stacklevel=0
        for k in range (len(a)):
            if (string.find(a[k][1],'ipython console') >0):
                stacklevel=k
        myf=sys._getframe(stacklevel).f_globals
        tb=myf['tb']
        return (tb)

# Get data if we do not have pyrap.tables using casa
def dget (vis, tel1, tel2):
    os.system('taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' % (vis, tel1, tel2, 'cl_temp.ms'))
    tb = define_tb()
    tb.open('cl_temp.ms')
    ut = tb.getcol('TIME')
    spw = tb.getcol('DATA_DESC_ID')
    uvw = tb.getcol('UVW')
    uvw = np.take(uvw,np.argwhere(spw==0),axis=1)[:,:,0]  # !!!!
    uvw = np.swapaxes(uvw,0,1)
    d = tb.getcol('DATA')
    tb.close()
    if spw.sum():    # (pol,chan,time+spw/cycl) -> (pol,chan,time)
        for i in np.unique(spw):
            new = np.take(d,np.argwhere(spw==i),axis=2)[:,:,:,0]
            try:
                d_out = np.concatenate((d_out,new),axis=1)
            except:
                d_out = np.copy(new)
        d = d_out
    return d,ut,uvw

