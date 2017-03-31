#  closure.py v1  Neal Jackson 2017 Jan 4
# Given a visibility file, return the scatter on the closure phases for a
# high signal-to-noise triangle. Scatter is 1.64 for random closure phases
# and about 1.00 for reasonably strong sources.
#
# I don't have pyrap.tables, and local IT cannot install it, so this 
# contains a CASA workaround this means that the output from pyrap.tables 
# is put into an ndarray in the same order that CASA tb.getcol gets it out. 
# I need this in order to do testing/debugging.

import sys,os,idx_tels,numpy as np, matplotlib, time
from matplotlib import pyplot as plt
from model_engine import *
try:
    import pyrap; from pyrap import tables as pt
    have_tables = True
except:
    have_tables = False

def closure (vis, tel, lastv=-1, plotfile='clplot.png'):
    # Find which number is which antenna
    itels = np.sort(idx_tels.get_idx_tels (vis, tel))
    if itels == []:
        return -1

    # Make three reference MSs with pointers
    print 'itels',itels
    if have_tables:
        d1,ut1,uvw = dget_t (vis,itels[0],itels[1])
        d2,ut2,uvw = dget_t (vis,itels[1],itels[2])
        d3,ut3,uvw = dget_t (vis,itels[0],itels[2])
    else:
        d1,ut1,uvw = dget_c (vis,itels[0],itels[1])
        d2,ut2,uvw = dget_c (vis,itels[1],itels[2])
        d3,ut3,uvw = dget_c (vis,itels[0],itels[2])
    a1,p1 = getap(d1[:lastv])
    a2,p2 = getap(d2[:lastv])
    a3,p3 = getap(d3[:lastv])
    clph = p1+p2-p3
    np.putmask(clph,clph>np.pi,clph-2*np.pi)
    np.putmask(clph,clph<-np.pi,clph+2*np.pi)
    # return a statistic - is 1.64 for random closure phase, less for coherent
    if len(plotfile):
        plt.plot(clph,'b+')
        plt.savefig(plotfile)
    return np.mean(np.gradient(np.unwrap(clph))**2)
