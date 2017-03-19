#--- model_engine.py: makes a (currently) two-component model
#    from u-v data on a particular closure triangle
#           v.1 Neal Jackson, 2015.09.29
#           v.2 NJ, 2017.01.16 converted to CASA, many changes
#           v.3 NJ, 2017.03.07 parallelised
import astropy,numpy as np,scipy,sys,os,glob,warnings,pyrap,multiprocessing,matplotlib
from pyrap import tables as pt; from matplotlib import pyplot as plt
from scipy import ndimage,optimize; from correlate import *; from mkgauss import *
plt.rcParams['image.origin'],plt.rcParams['image.interpolation']='lower','nearest'
warnings.simplefilter('ignore')
RAD2ARC,LIGHT,FIRSTNPY = 3600.*180./np.pi, 2.99792458E+8, './first_2008.simple.npy'
MX,MY,MF,MW,MR,MP = range(6)

# Make a movie from a set of png files
def movie (fps=4):
    a=np.sort(glob.glob('model_engine_*.png'))
    for i in range(len(a)):
        os.system('convert %s model_engine%03d.jpeg'%(a[i],i))
    command = "mencoder \"mf://model_engine*.jpeg\" -mf fps=%d -o model_engine.avi -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=4800" %fps
    os.system(command)
    os.system('rm model_engine*.jpeg')

# Work out real and imaginary parts of a list of visibilities given model and uvw's
def uvw2reim (uvw, model):
    u,v,re,im = uvw[:,0],uvw[:,1],0.,0.
    for m in model:     # note - sign in m[MX] as this is HA not RA
        cmpamp,cmpphs = m[MF],2.*np.pi*(-u*m[MX]+v*m[MY])/RAD2ARC
        if m[MW]!=0.0:
            sphi,cphi = np.sin(np.deg2rad(m[MP])),np.cos(np.deg2rad(m[MP]))
            tc = np.pi*(m[MW]/RAD2ARC)*np.hypot(v*cphi+u*sphi,m[MR]*(u*cphi-v*sphi))
            cmpamp = m[MF]*np.exp(-0.3696737602*tc*tc)
        re += cmpamp * np.cos(cmpphs)
        im += cmpamp * np.sin(cmpphs)
    return re,im

# From a MS, return the index numbers of telescopes given as an array of 3 strings
def get_idx_tels (data, tel):
    command = 'taql \'select NAME from %s/ANTENNA\' >closure_which'%data
    os.system(command)
    f = open('closure_which')
    for line in f:
        if not 'select' in line:
            try:
                a = int(line) # it's just a number, use STATION column instead
                f.close()
                command = 'taql \'select STATION from %s/ANTENNA\' >closure_which'%data
                os.system(command)
                break
            except:
                f.close()
                break
    idx_tels, iline = [-1,-1,-1], 0
    f = open('closure_which')
    for line in f:
        if not 'select' in line:
            for i in range(3):
                if tel[i]==line[:len(tel[i])]:
                    idx_tels[i] = iline
            iline += 1
    f.close()
    if -1 in idx_tels:
        print 'Did not find one or more of the telescopes'
        return []
    return idx_tels

# add up all the real and imaginary for all channels and spw, assuming d[time,spw]
def squash_amp_clph (d,spw,nspw,pol=0):
    a = np.zeros((len(d)/nspw),dtype='complex')
    n = np.zeros((len(d)/nspw))
    for i in range(len(d)):
        idx = i/nspw
        a[idx] += np.nansum(d[i]['DATA'][:,pol])
        n[idx] += len(~np.isnan(d[0]['DATA'][:,pol]))
    return abs(a/n),norm(np.arctan2((a/n).imag,(a/n).real),False,False)

# Return two amplitude arrays and a closure phase array from 3 data arrays
def get_amp_clph (d01,d02,d12,spw,otel,nspw=0,pol=0):
    global cp012,a01,a02,a12
    sp = np.array([])
    for i in range(len(d01)):
        sp = np.append(sp, spw[i]['DATA_DESC_ID'])
    nspw = len(np.unique(sp))
    a01,p01 = squash_amp_clph (d01,spw, nspw)
    a02,p02 = squash_amp_clph (d02,spw, nspw)
    a12,p12 = squash_amp_clph (d12,spw, nspw)
    cp012 = np.rad2deg(p01*otel[0]-p02*otel[1]+p12*otel[2])
    np.putmask(cp012,cp012>180.0,cp012-360.0)
    np.putmask(cp012,cp012<-180.0,cp012+360.0)

def get_uvw_table (t):
    for u in t.select('UVW'):
        try:
            uvw = np.vstack ((uvw,u['UVW']))
        except:
            uvw = np.copy(u['UVW'])
    return uvw

# Get data and u-v arrays from a measurement set on a given triangle
def data_extract (data):
    global uvw01,uvw02,uvw12
    ctemplate = 'taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving cl_tmp%d.ms\''
    os.system ('msoverview in='+data+' >temp_msoverview;cat temp_msoverview')
    f = open('temp_msoverview')
    nif = 0
    for l in f:    # parse msoverview output to get channel information
        print l
        if 'TOPO' in l:
            nif+=1
            nchan = int(l.split()[2])
            ch0 = float(l.split()[4])
            chw = float(l.split()[5])
        if 'J2000' in l:    # and coordinate information
            sra,sdec = l[32:].split()[0].split(':'),l[32:].split()[1].split('.')
            isneg = -1.0 if sdec[0][0]=='-' else 1.0
            sra =  np.asarray(sra,dtype='float')
            sdec = np.asarray(sdec,dtype='float')
            ra = 15.*sra[0]+sra[1]/4.+sra[2]/240.
            dec = sdec[0]+sdec[1]/60.+sdec[2]/3600.
    f.close(); os.system ('rm temp_msoverview')
    wlength = np.mean(LIGHT/(ch0*1.e6 + chw*1.e3*np.arange(nif*nchan)))
    itel = np.array(get_idx_tels (data,trname))
    # order of telescopes on baselines 0-1, 0-2, 1-2; -1 if in "wrong" order
    otel = 1.-2.*np.asarray([itel[0]>itel[1],itel[0]>itel[2],itel[1]>itel[2]],dtype=float)
    btel = [min(itel[0],itel[1]),max(itel[0],itel[1]),min(itel[0],itel[2]),\
            max(itel[0],itel[2]),min(itel[1],itel[2]),max(itel[1],itel[2])]
    os.system (ctemplate % (data,btel[0],btel[1],1))
    os.system (ctemplate % (data,btel[2],btel[3],2))
    os.system (ctemplate % (data,btel[4],btel[5],3))
    t01,t02,t12 = pt.table('cl_tmp1.ms'),pt.table('cl_tmp2.ms'),pt.table('cl_tmp3.ms')
    spw = t01.select('DATA_DESC_ID')
    d01,d02,d12 = t01.select('DATA'), t02.select('DATA'), t12.select('DATA')
    get_amp_clph(d01,d02,d12,spw,otel)
    times = np.array([])
    for t in t01.select('TIME'):
        times = np.append(times,t['TIME'])
    uvw01 = get_uvw_table(t01)/wlength      # km->wavelength
    uvw02 = get_uvw_table(t02)/wlength
    uvw12 = get_uvw_table(t12)/wlength
    print trname,'-> antenna numbers:',itel
    print 'Baseline lengths: %s-%s: %dkm %s-%s: %dkm %s-%s: %dkm' % \
       (trname[0],trname[1],int(np.sqrt((uvw01[0]**2).sum())*wlength/1000),\
        trname[0],trname[2],int(np.sqrt((uvw02[0]**2).sum())*wlength/1000),\
        trname[1],trname[2],int(np.sqrt((uvw12[0]**2).sum())*wlength/1000))
    os.system('rm -fr cl_tmp*.ms')
    return itel,np.mean(wlength),ra,dec

# --------------------------------------------------

def model_extract (model,itel):
    otel = 1.-2.*np.asarray([itel[0]>itel[1],itel[0]>itel[2],itel[1]>itel[2]],dtype=float)
    re01,im01 = uvw2reim (uvw01,model)
    re12,im12 = uvw2reim (uvw12,model)
    re02,im02 = uvw2reim (uvw02,model)
    ph01 = norm(np.rad2deg(np.arctan2(im01,re01)))
    ph02 = norm(np.rad2deg(np.arctan2(im02,re02)))
    ph12 = norm(np.rad2deg(np.arctan2(im12,re12)))
    clph = norm(ph01*otel[0] - ph02*otel[1] + ph12*otel[2])
    return np.hypot(re01,im01),np.hypot(re02,im02),clph

def plotimg (A01,A02,CP012,model,goodness,itel,aplot,gcou):
    ells = []
    for i in model:
        if i[MW]!=0.0:
            ells.append(matplotlib.patches.Ellipse(xy=[i[MX],i[MY]],width=2.*i[MW],\
                height=2.*i[MW]*i[MR],lw=0.5,angle=i[MP]+90.,fill=0.0))
        else:
            ells.append(matplotlib.patches.Ellipse(xy=[i[MX],i[MY]],width=0.1,\
                height=0.1,lw=0.5,angle=0.,fill=1.0,color='red'))
    if plottype in [1,10]:
        plt.subplot2grid((6,5),(0,0),rowspan=2)
        scalarray = np.sort(np.ravel(a01)); ls = len(scalarray)
        amp_min, amp_max = scalarray[int(0.1*ls)], scalarray[int(0.9*ls)]
        plt.imshow(a01,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(0,1),rowspan=2)
        plt.imshow(A01,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(0,2),rowspan=2)
        plt.imshow(a01-A01,aspect='auto',\
                   vmin=2*amp_min-amp_max,vmax=2*amp_max-amp_min)
        plt.subplot2grid((6,5),(2,0),rowspan=2)
        scalarray = np.sort(np.ravel(a02)); ls = len(scalarray)
        amp_min, amp_max = scalarray[int(0.1*ls)], scalarray[int(0.9*ls)]
        plt.imshow(a02,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(2,1),rowspan=2)
        plt.imshow(A02,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(2,2),rowspan=2)
        plt.imshow(a02-A02,aspect='auto',\
                   vmin=2*amp_min-amp_max,vmax=2*amp_max-amp_min)
        plt.subplot2grid((6,5),(4,0),rowspan=2,yticks=[])
        plt.imshow(cp012,aspect='auto',vmin=-180,vmax=180)
        plt.subplot2grid((6,5),(4,1),rowspan=2,yticks=[])
        plt.imshow(CP012,aspect='auto',vmin=-180,vmax=180)
        plt.subplot2grid((6,5),(4,2),rowspan=2,yticks=[])
        plt.imshow(CP012-cp012,aspect='auto')
        plt.subplot2grid((6,5),(3,3),rowspan=3,colspan=3,xticks=[],yticks=[])
        plt.imshow(aplot,cmap=matplotlib.cm.gray_r)
        plt.title('X2=%.3f'%goodness)
    elif plottype in [2,20]:
        plt.subplot2grid((6,5),(0,0),rowspan=2,colspan=3)
        plt.plot(a01,'b-'); plt.plot(A01,'r-')
        plt.legend([trname[0]+'-'+trname[1]],fontsize=6)
        plt.subplot2grid((6,5),(2,0),rowspan=2,colspan=3)
        plt.plot(a02,'b-'); plt.plot(A02,'r-')
        plt.legend([trname[0]+'-'+trname[2]],fontsize=6)
        plt.subplot2grid((6,5),(4,0),rowspan=2,colspan=3)
        plt.plot(cp012,'b-'); plt.plot(CP012,'r-')
        plt.subplot2grid((6,5),(4,3),rowspan=2,colspan=3,xticks=[],yticks=[])
        plt.imshow(aplot,cmap=matplotlib.cm.gray_r)
        plt.title('X2=%.3f'%goodness)
    ax = plt.subplot2grid((6,5),(0,3),rowspan=3,colspan=2)
    for i in range(len(ells)):
        e = ells[i]
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
    ax.set_xlim(-glim,glim);ax.set_ylim(-glim,glim)
    plt.grid()
    if gcou==-1:   # find how many png files, make the next one
        gcou = len(glob.glob('model_engine_*.png'))
    plt.title('model_engine_%03d.png'%gcou)
    if plottype in [10,20]:
        plt.savefig('model_engine_%03d.png'%gcou)

#===================================================
def ndiff (a,b):
    sqd = np.array([])
    for i in range(-len(a)/2,len(a)/2):
        a1 = np.roll(a,i)
        idx1,idx2 = max(0,i),min(len(a),len(a)+i)
        sqd = np.append(sqd, np.mean(((b-a1)[idx1:idx2]**2)))
    return sqd

# this is the bit that may need replacing - maybe fit splines and look at the spline points?

def get_goodness(A01,A02,CP012):
    beta = 0.00001    #   this is a pretty vital parameter
    ascat = np.median(abs(np.gradient(np.ravel(a02))))    
    cscat = np.median(abs(np.gradient(np.ravel(cp012))))
    sq = ndiff(a01*np.nanmean(A01)/np.nanmean(a01),A01)/ascat**2 + \
         ndiff(a02*np.nanmean(A02)/np.nanmean(a02),A02)/ascat**2 + \
         ndiff(cp012,CP012)/cscat**2
    difmin = 0.5*len(a01)-np.argwhere(sq==sq.min())[0][0]
    return sq.min() + beta*difmin**2

def mod_func (x0, *x):
    model,opt,itel,aplot,gcou,iy,ix = x
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = x0
    model = model.reshape(len(model)/6,6)
    A01,A02,CP012 = model_extract (model,itel)
    if ampfiddle:
        A01 *= np.median(a01)/np.median(A01)
        A02 *= np.median(a02)/np.median(A02)
    goodness = get_goodness(A01,A02,CP012)
    if plottype:
        plotimg (A01,A02,CP012,model,goodness,itel,aplot,gcou)
        if plottype in [1,2]:
            plt.draw()
            plt.pause(0.001)
            plt.clf()
    return goodness

def getmodel(coord,beam):
    gsiz = int(gridsize/(bsub*beam))   # original grid very big
    ginc = bsub*beam                   # arcsec/pix grid size
    grid = np.ones((gsiz,gsiz))*np.nan # original grid full of NaN
    first = np.load(FIRSTNPY)
    pflux,pcoord = [],[]
    a = correlate (np.array([coord]),0,1,first,0,1,0.5/60)
    cosdec = np.cos(np.deg2rad(coord[1]))
    print 'Found: %d FIRST sources'%len(a)
    for i in range(len(a)):
        fi = first[int(a[i,1])]
        scoord = astropy.coordinates.SkyCoord(fi[0],fi[1],unit='degree')
        scoord = scoord.to_string(style='hmsdms')
        scoord = scoord.replace('d',':').replace('m',':')
        scoord = scoord.replace('h',':').replace('s','')
        print '%s flux=%.1f shape=(%.1fx%.1f PA%.1f)'%(scoord,fi[2],fi[4],fi[5],fi[6])
        pix = 0.5*gsiz+3600.*(fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        try:
            pcoord.append(fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        except:
            pcoord = (fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        pflux.append(fi[2])
        new = mkgauss([gsiz,gsiz],pix,fi[2],fi[4]/ginc,fi[5]/fi[4],fi[6])
        np.putmask(grid,new>0.004*new.max(),i)  # decrease 0.01 if missing cpts
    # grid now consists of original very big grid, with numbers instead of NaN where
    # the secondary might be
    if pflux:    # shrink the grid around the Gaussian near FIRST sources
        while all(np.isnan(grid[0,:])) and all(np.isnan(grid[-1,:])) and \
              all(np.isnan(grid[:,0])) and all(np.isnan(grid[:,-1])):
            grid = grid[1:-1,1:-1]
    return grid,ginc,grid.shape[0],pflux,pcoord

def parallel_function(f):
    def easy_parallize(f, sequence):
        from multiprocessing import Pool
        ncores = max(1,int(os.popen('nproc').read())-1)  # use all cores - 1
        pool = Pool(processes=ncores) # depends on available cores
        result = pool.map(f, sequence) # for i in sequence: result[i] = f(i)
        cleaned = [x for x in result if not x is None] # getting results
        cleaned = np.asarray(cleaned)
        pool.close() # not optimal! but easy
        pool.join()
        return cleaned
    from functools import partial
    return partial(easy_parallize, f)

def grid_search_thread (k):
    a =  mod_func(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7])
    f=open('grid_search_thread.log','a')
    f.write('%d %d %f\n'%(k[-2],k[-1],a))
    f.close()
    return a

def grid_search (model,cpt,gridcpt,itel,aplot,gcou,grid,gsiz,ginc):
    opt = np.zeros_like(model,dtype='bool')
    args = []
    import copy
    for ix in range(gsiz):
        x = ginc*(ix-gsiz/2.0)   # x,y in arcsec; a in ginc-size pixels
        for iy in range(gsiz):
            if grid[iy,ix] == gridcpt:
                y = ginc*(iy-gsiz/2.0)
                model[cpt][0],model[cpt][1] = x,y
                arg = ([],model,opt,itel,aplot,gcou,iy,ix)
                args.append(copy.deepcopy(arg))   # otherwise overwrites all elements
                gcou += 1
    print 'Starting grid search with',len(args),'points'
    os.system('rm grid_search_thread.log')
    grid_search_thread.parallel = parallel_function(grid_search_thread)
    parallel_result = grid_search_thread.parallel (args)
    for i in range(len(parallel_result)):
        aplot[args[i][-2],args[i][-1]] = parallel_result[i]
    np.putmask(aplot,np.isnan(aplot),np.nanmax(aplot))
    model[cpt,:2] = ginc*(np.asarray(ndimage.measurements.minimum_position \
            (aplot)[::-1])-0.5*np.asarray(grid.shape))
    return model,aplot

def recentroid (model,startcpt,endcpt,ginc,gsiz):
    model[startcpt:endcpt+1,:2]/=ginc
    model[startcpt:endcpt+1,3]/=ginc
    centroid = np.array([0.0,0.0])
    for i in range(startcpt,endcpt+1):
        centroid += np.array(model[i,2]*model[i,:2])
    centroid /= np.sum(model[startcpt:endcpt+1,2])
    model[startcpt:endcpt+1,:2] += np.array([0.5*gsiz,0.5*gsiz]) - centroid
    return model

def adjust_all (model,itel,aplot,gcou):
    opt = np.ones_like(model,dtype='bool')
    x0 = np.ravel(model)[np.ravel(opt)]
    args = (model,opt,itel,aplot,-1,0,0)
    xopt = optimize.fmin(mod_func, x0, args=args, maxiter=100)
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = xopt
    model = model.reshape(len(model)/6,6)
    return model

def refine_points (model,cpts,itel,aplot,gcou):
    opt = np.zeros_like(model,dtype='bool')
    for i in cpts:
        opt[i,2:] = True
    x0 = np.ravel(model)[np.ravel(opt)]
    args = (model,opt,itel,aplot,-1,0,0)
    xopt = optimize.fmin(mod_func, x0, args=args, maxiter=20)
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = xopt
    model = model.reshape(len(model)/6,6)
    return model

def norm(a,isdeg=True,isred=True):
    nlim = 180.0 if isdeg else np.pi
    a = a%(2*nlim) if isred else a
    np.putmask(a,a>nlim,a-2.*nlim)
    np.putmask(a,a<-nlim,a+2.*nlim)
    return a

def write_skymodel (ra,dec,model,outname):
    if outname!='':
        f = open(outname,'w')
        f.write ('# (Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, Orientation) = format\n')
    for i in range(len(model)):
        cosd = 3600.*np.cos(np.deg2rad(dec))
        s = astropy.coordinates.SkyCoord(ra-model[i,0]/cosd,dec+model[i,1]/3600,unit='degree')
        s = s.to_string(style='hmsdms')
        sra = s.split()[0].replace('h',':').replace('m',':').replace('s','')
        sdec = s.split()[1].replace('d','.').replace('m','.').replace('s','')
        print 'ME%d, GAUSSIAN, %s, %s, %f, %f, %f, %f'%(i,sra,sdec,model[i,MF],\
              model[i,MW],model[i,MW]*model[i,MR],np.rad2deg(model[i,MP]))
        if outname!='':
            f.write('ME%d, GAUSSIAN, %s, %s, %f, %f, %f, %f\n'%(i,sra,sdec,model[i,MF],\
                  model[i,MW],model[i,MW]*model[i,MR],np.rad2deg(model[i,MP])))
    if outname!='':
        f.close()

#============= main script ========================
def mainscript(data,TRNAME,BSUB=0.3,GRIDSIZE=12.0,PLOTTYPE=20,AMPFIDDLE=True,outname='model_engine.sky'):
    global bsub,gridsize,plottype,ampfiddle,glim,trname
    bsub,gridsize,plottype,ampfiddle,trname,gcou = BSUB,GRIDSIZE,PLOTTYPE,AMPFIDDLE,TRNAME,0
    os.system('rm model_engine*.png')
    itel,wv,ra,dec = data_extract (data)
    s_amp = np.sort(np.ravel(a01)); ls = len(s_amp)
    flux1,flux2 = np.median(s_amp), 0.5*(s_amp[int(0.99*ls)]-s_amp[int(0.01*ls)])
    flux = np.median(a01)
    beam = RAD2ARC/np.nanmax(abs(uvw01)); print 'Beam:',beam,'arcsec'
    grid,ginc,gsiz,pflux,pcoord = getmodel (np.array([ra,dec]),beam)
    if not len(pflux):   # no FIRST source, search the whole grid
        grid = np.zeros_like (grid)
    aplot = np.ones_like(grid)*np.nan
    glim = 0.5*ginc*gsiz
    if plottype in [1,2]:
        plt.ion()
        fig = plt.figure()
    if len(pflux) < 2:    # zero, or one FIRST source
        model = np.array([[0.0,0.0,flux1,0.5,1.0,0.0],[-1.0,-2.0,flux2,0.5,1.0,0.0]])
        model,aplot = grid_search (model,1,0,itel,aplot,gcou,grid,gsiz,ginc)   
# fix one cpt, grid-search posn of 2nd.
        print 'Model after grid search:'
        write_skymodel (ra,dec,model,'')
        model[0,MW] = model[1,MW] = 1.0
        model = refine_points (model,[0,1],itel,aplot,gcou)  
# fit for size/orientation only
        print 'Model after refining points:'
        write_skymodel (ra,dec,model,'')
        model = recentroid (model,0,1,ginc,gsiz)  # put centroid in centre of image
        print 'Model after recentroiding:'
    else:   # >1 FIRST source - not tested yet
        model = np.zeros((len(pflux),6))
        for i in range(len(pflux)):
            model[i,:2] = pcoord[i]
            model[i,3] = pflux[i]*flux/pflux.sum()
            model[i,4:] = [0.5,1.0,0.0]
        model = adjust_all(model,itel,aplot,gcou)
    if plottype in [10,20]:
        movie()
    write_skymodel (ra,dec,model,'')
    write_skymodel (ra,dec,model,outname)
#mainscript('./SIM5.ms', ['ST001','DE601','DE605HBA'])
mainscript('./PP1_av_corr.ms', ['ST001','DE601HBA','DE605HBA'])
