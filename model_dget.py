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

def taql_get (vis, ant1, ant2, fname):
    print 'Selecting ants',ant1,ant2
    os.system('taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' % (vis, ant1, ant2, fname))

# Turn array in (pol,chan,time+spw cycle) into (pol,chan,time) in fq order
def spec_flatten (d, spw, isuv=False):
    for i in np.unique(spw):
        new = np.take(d,np.argwhere(spw==i),axis=2)[:,:,:,0]
        try:
            d_out = np.concatenate((d_out,new),axis=1)
        except:
            d_out = np.copy(new)
    return d_out

# Get data if we do not have pyrap.tables using casa
def dget (vis, tel1, tel2):
    taql_get (vis, tel1, tel2, 'cl_temp.ms')
    tb = define_tb()
    tb.open('cl_temp.ms')
    ut = tb.getcol('TIME')
    spw = tb.getcol('DATA_DESC_ID')
    uvw = tb.getcol('UVW')
    uvw = np.take(uvw,np.argwhere(spw==0),axis=1)[:,:,0]  # !!!!
    uvw = np.swapaxes(uvw,0,1)
    d = tb.getcol('DATA')
    tb.close()
    if spw.sum():
        d = spec_flatten (d,spw)
    return d,ut,uvw

