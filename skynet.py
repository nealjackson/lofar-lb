import numpy as np,os,sys,glob,astropy
import closure,model_engine

# Requires:
# 1) list of MSs of sources within the field, with delay calibration
# 2) list of coordinates of those fields (np.array([[ra,dec],[ra,dec]....])
# 3) list of three closure telescopes (e.g. ['TS001','DE601','DE605'])

# 4ch/8s 20Gb per source (FOV) 1ch/16s 3 GB/source to go to 5' fields
# eor scripts gives out cal table as a lofar parmdb
# use ndppp to apply parmdb to MS

CTHR = 1.3

def ndppp_vs_model (vis, mod, solint):
    f=open('NDPPP.parset')
    f.write('msin='+vis+'\n')
    f.write('msin.datacolumn=DATA\n')
    f.write('steps=[gaincal]\n')
    f.write('gaincal.parmdb='+vis+'/instrument\n')
    f.write('gaincal.caltype=phaseonly\n')
    f.write('gaincal.solint='+solint+'\n')
    f.write('gaincal.maxiter=100\n')
    f.close()
    os.system('NDPPP NDPPP.parset')

def sky_eval (flist, coordlist, closure_tels):
    ns = len(flist)
    if ns != len(coordlist):
        print 'File list must be same length as coord list'
        return
    if not ns:
        print 'No sources in field'
        return
    closure_stats = np.array([])
    for i in range(ns):
        closure_this = closure.closure(flist[i],closure_tels)
        closure_stats = np.append(closure_stats, closure_this)
    # make models and solutions for all closure_stats < 1.3
        if closure_this < CTHR:
            model_engine.mainscript (flist[i],plottype=0,skymodel=flist[i]+'_mod')
            sanity_check_model (flist[i]+'_mod')
            ndppp_vs_model (flist[i],flist[i]+'_mod',15)
            eormap (flist[i],parmdb=flist[i]+'_parmdb')
            sanity_check_corrections (flist[i]+'_parmdb')
    
    
    
    
