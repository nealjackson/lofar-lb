#!/usr/bin/env python
import sys,os,time,numpy as np,pyrap,matplotlib
from pyrap import tables as pt
from matplotlib import pyplot as plt

#  closure.py v5  Neal Jackson 2018 May 28
# Given a set of visibility files and closure triangles, return the significance of 
# detection of each source on each triangle as an array. Optionally make .npy dumps
# of each of the closure phase dynamic spectra. See LBWG memo number 3 for details.
#
# arguments:
#    vis = list of visibility file names
#    tel = np.array nx3 strings of telescopes to use
#    lastv = only process this number of data points (default all)
#    bchan, echan: first and last channel to use (NB channels numbered from
#       channel 0 of first spw to last channel of last spw, ie in freq order)
#    pol = polarization (default first in file, RR)
#    fromaips: if file was written by AIPS, will have tel names in STATIONS column
#
# returns an mxn array of scatters, where n is the number of closure triangles and m the number of MSs
# These are the significance in sigma of the detection

def closure(vis,tel,lastv=-1,pol=0,bchan=0,echan=-1,fromaips=False,dump=False):

    # Find target source id

    target_id = vis.split('/')[-1].split('_')[0]
    print "target_id", target_id
    print "Antennas for closure phase", tel

    # Find array of requested telescopes and list of telescopes in data

    antcol = 'STATION' if fromaips else 'NAME'
    command = 'taql \'select %s from %s/ANTENNA\' >closure_txt'%(antcol,vis)
    os.system(command)
    os.system('grep -v select closure_txt >closure_which')
    idxtel = np.loadtxt('closure_which',dtype='S')
    atel = np.unique(np.ravel(tel))

    # For each requested telescope, determine its position in the list
    # If more than one telescope in the list match to within the number
    #   of letters in the requested telescope, keep the first (so CS002
    #   will match to CS002HBA0 and CS002HBA1 will be ignored)
    # Keep a list of telescopes not found, to print if we need to crash

    notfound = []
    aidx = np.array([],dtype='int')
    for a in atel:
        found_this = False
        for i in range(len(idxtel)):
            if a==idxtel[i][:len(a)]:
                aidx = np.append(aidx,i)
                found_this = True
        if not found_this:
            notfound.append (a)
    if len(notfound):
        print 'The following telescopes were not found:',notfound
    aidx_s = np.sort(aidx)

    # Make a smaller MS 'as plain' with the required baselines. This is slow
    # but only needs doing once for an arbitrary number of baselines.

    if os.path.exists ('cl_temp.ms'):
        os.system('rm -fr cl_temp.ms')
    command = 'taql \'select from %s where ' % vis
    for i in range (len(aidx_s)):
        for j in range (i+1, len(aidx_s)):
            command += ('ANTENNA1==%d and ANTENNA2==%d' % \
                            (aidx_s[i],aidx_s[j]))
            if i==len(aidx_s)-2 and j==len(aidx_s)-1:
                command += (' giving cl_temp.ms as plain\'')
            else:
                command += (' or ')
    print 'Selecting smaller MS cl_temp.ms, this will take about 4s/Gb:'
    print command
    os.system (command)

    # Loop around the requested closure triangles

    clstats = np.array([])
    for tr in tel:
        print 'Doing triangle:',tr
        tri = np.array([],dtype='int')
        for i in range(3):
            tri = np.append (tri, aidx[np.argwhere(atel==tr[i])[0][0]])
        tri = np.sort(tri)

        # Make three reference MSs with pointers into the small MS

        command = 'taql \'select from cl_temp.ms where ANTENNA1==%d and ANTENNA2==%d giving closure_temp1.ms\'' %(tri[0],tri[1])
        os.system(command)
        command = 'taql \'select from cl_temp.ms where ANTENNA1==%d and ANTENNA2==%d giving closure_temp2.ms\'' %(tri[1],tri[2])
        os.system(command)
        command = 'taql \'select from cl_temp.ms where ANTENNA1==%d and ANTENNA2==%d giving closure_temp3.ms\'' %(tri[0],tri[2])
        os.system(command)

        # Load data arrays and get amp, closure phase

        dtab = [pt.table('closure_temp1.ms').select('DATA'),\
                pt.table('closure_temp2.ms').select('DATA'),\
                pt.table('closure_temp3.ms').select('DATA')]
        ttab = [pt.table('closure_temp1.ms').select('TIME'),\
                pt.table('closure_temp2.ms').select('TIME'),\
                pt.table('closure_temp3.ms').select('TIME')]
        stab = [pt.table('closure_temp1.ms').select('DATA_DESC_ID'),\
                pt.table('closure_temp2.ms').select('DATA_DESC_ID'),\
                pt.table('closure_temp3.ms').select('DATA_DESC_ID')]        
        d = [np.array([]),np.array([]),np.array([])]
        t = [np.array([]),np.array([]),np.array([])]
        s = [np.array([]),np.array([]),np.array([])]
        for i in range(3):      # assumes file regular!
# If crashes here, need to check that all SB/channels are present for each BL
            for j in range(len(dtab[i])):
                new = dtab[i][j]['DATA'][:,pol]
                try:
                    d[i] = np.vstack((d[i],new))
                except:
                    d[i] = np.copy(new)
            for j in range(len(ttab[i])):
                t[i] = np.append(t[i],ttab[i][j]['TIME'])
            for j in range(len(stab[i])):
                s[i] = np.append(s[i],stab[i][j]['DATA_DESC_ID'])
            nspw = len(np.unique(s[i]))
            time_ints = d[i].shape[0]/nspw
            d[i] = d[i].reshape(time_ints,nspw*d[i].shape[1])
            d[i] = np.arctan2(d[i].imag,d[i].real)
            np.putmask(d[i],d[i]==0.0,np.nan)

        from scipy.fftpack import *
        p = np.asarray(d[0]+d[1]-d[2],dtype='float')
        p = p[:,bchan:p.shape[1] if echan == -1 else echan]
        np.putmask(p,p>np.pi,p-2.*np.pi)   # no need to worry about nans
        np.putmask(p,p<-np.pi,p+2.*np.pi)
        if dump:
            np.save('clph_%s_%s_%s_%s' % (vis,tr[0],tr[1],tr[2]),p)
        prand = np.random.random(p.shape[0]*p.shape[1]).reshape(p.shape)*2*np.pi-np.pi
        np.putmask(p,np.isnan(p),prand)
        pp = np.cos(p)+1j*np.sin(p)
        pprand = np.cos(prand)+1j*np.sin(prand)
        qbig = fftshift(fft2(pp))
        qbigrand = fftshift(fft2(pprand))
        FFSIZE = min(30,d[0].shape[1]/2)
        yr = (qbig.shape[0]/2 - FFSIZE,qbig.shape[0]/2+FFSIZE)
        xr = (qbig.shape[1]/2 - FFSIZE,qbig.shape[1]/2+FFSIZE)
        q = qbig[yr[0]:yr[1],xr[0]:xr[1]]
        qrand = qbigrand[yr[0]:yr[1],xr[0]:xr[1]]
        qq = np.sort(np.ravel(abs(q)))
        qqrand = np.sort(np.ravel(abs(qrand)))
        x,y = np.arange(-FFSIZE,FFSIZE),np.arange(-FFSIZE,FFSIZE)
        xx,yy = np.meshgrid(x,y)
        z = np.hypot(xx,yy)
        clthis = 0.0
        for i in range(3,FFSIZE):
            zz=np.zeros_like(z)
            np.putmask(zz,z<i,1.0)
            sig = ((zz*abs(q)).sum()-(zz*abs(qrand)).sum()) \
                      /(np.std(abs(qrand))*np.sqrt(zz.sum()))
            clthis = max(clthis,sig)
        plt.subplot(131,xticks=[],yticks=[]);plt.imshow(p,origin='lower',aspect=0.15,cmap=matplotlib.cm.jet)
        plt.subplot(132);plt.imshow(abs(qrand),origin='lower',vmax=4000,aspect='equal')
        plt.subplot(133);plt.imshow(abs(q),origin='lower',aspect='equal',vmax=4000);plt.show()
        os.system('rm -fr closure_temp*ms')
        if os.path.exists ('closure_which'):
           os.system('rm closure_which')
        clstats = np.append (clstats, clthis)
    return clstats


def main(ms_input,station_input,lastv=-1,pol=0,bchan=0,echan=-1,fromaips=False,dump=False):
    """
    Deriving closure phases of all directions
   
    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile
    station_input : List of stations for closure triangles, separated by
        semicolons. There must be nx3 of those.
    lastv (int):   Index of last visibility
    pol     (int, default 0): Polarization number
    bchan   (int, default 0): First channel number
    echan   (int, default -1): Last channel number
    fromaips If True, file converted from aips (uses STATION for ant. names)
    dump     If True,dump a lot of .npy files of closure phase dynamic spectra
        
    Returns
    -------
    result : dict
        Dict with the name of the generated mapfile
    """
    
    mslist = str(ms_input).lstrip('[').rstrip(']').replace("'","").replace(" ","").split(',')
    stationlist = str(station_input).lstrip('[').rstrip(']').replace("'","").replace(" ","").split(';')
    tel = np.asarray(stationlist)
    if len(tel)%3:
        print 'Station list is not a multiple of 3 long, truncating it...'
        tel = tel[0:3*(len(tel)//3)]
    tel = tel.reshape(len(tel)/3,3)
    for ms in mslist:
        print 'Now operating on', ms
        direction = ms.split('/')[-1].split('_')[0]
        scatter_cp = closure(ms,tel,lastv=lastv,pol=pol,\
                             bchan=bchan,echan=echan,fromaips=fromaips,dump=dump)
        infostring = '\n Scatter for the direction %s is %s \n' % (direction,scatter_cp[0])
        print infostring
        command = 'echo \'%s\' >> closure_phases.txt' % infostring
        os.system(command)
        try:
            allcp = np.vstack((allcp,scatter_cp))
        except:
            allcp = np.copy(scatter_cp)
    return allcp

