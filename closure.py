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

def closure (vis='L519076.ms', tel=['ST001','DE601','DE605'], lastv=-1, plotfile='clplot.png'):
    # Find which number is which antenna
    itels = np.sort(idx_tels.get_idx_tels (vis, tel))
    if itels == []:
        print 'Failed to find, returning'
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
    os.system('rm -fr closure_temp*ms')
    os.system('rm closure_which')
    if len(plotfile):
        plt.plot(clph,'b+')
        plt.savefig(plotfile)
    return np.mean(np.gradient(np.unwrap(clph))**2)


####--------------- I am working on the stuff below this line ----------------

# Given a list of files processed by the generic-pipeline into subbands, calculate
# delays using CASA gaincal (stub - will eventually be replaced either by CASA fringe
# fitter or transferring delay solutions from the calibrator)

#TESTMSLIST = ['/data/scratch/lb_bw/lbgp/working_directory/lbgp/L401323_SBgr000-10_PP1_uv.dppp.ndppp_prep_target','/data/scratch/lb_bw/lbgp/working_directory/lbgp/L401323_SBgr000-10_PP2_uv.dppp.ndppp_prep_target']

#def fringefit (mslist=TESTMSLIST):
#    closure_stats = []
#    for msthis in mslist:
#        closure_stats.append (closure(vis=msthis, tel=['DE601','DE605','ST001'], lastv=300))
#    print 'Best source appears to be'

# Test datasets: closure('L519076.ms',['ST001','DE601HBA','DE605HBA'])
#               closure('PP1_av_corr.ms',['ST001','DE601HBA','DE605HBA'])

