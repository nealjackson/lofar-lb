#  Interactive model-fitter  NJ 2017.08.23
#  Change the last two lines of the program and run
#
#  In the plotting surface:
#  Hit '0' to move the initial component around, '1' etc for new components
#  Use f,F,w,W,r,R,p,P to adjust the flux, width, axis ratio, pos ang of any
#     component (the one nearest to the click will be adjusted)
#  Aim is to match the amplitudes and closure phases plotted on the RHS
#  Output is stored as a skymodel file iengine.skymodel
#
import numpy as np,astropy,matplotlib,glob,os,sys,time,warnings,idx_tels
from matplotlib import pyplot as plt, patches as patches
from matplotlib.patches import Ellipse
from matplotlib.widgets import Slider, Button, RadioButtons  #, TextBox
from taql_funcs import *
from astropy import coordinates
try:
    import pyrap; from pyrap import tables as pt
    have_tables = True
except:
    have_tables = False
CASAPY = '/pkg/casa-release-4.7.0-1-el6/bin/casa'
#  Change the following things if needed:
PLOTSIZE = (28,22)    # size of plot on screen (32.5x25 is a big screen)
CSCALE = matplotlib.cm.hot
MX,MY,MF,MW,MR,MP = range(6)
warnings.simplefilter('ignore')
RAD2ARC,LIGHT,FIRSTNPY = 3600.*180./np.pi, 2.99792458E+8, './first_2008.simple.npy'
model = np.array([[0.0,0.0,1.0,0.4,0.0,0.0]])
trname = ['ST001','DE601','DE605']
plt.rcParams['image.origin']='lower'
plt.rcParams['image.interpolation']='nearest'
# ============= stuff from model_engine          ==============

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

def dget_c (vis, tel1, tel2):
    f,fo = open('model_dget.py'), open('temp.py','w')
    for line in f:
        fo.write(line)
    f.close()
    os.system('rm d.npy; rm ut.npy; rm uvw.npy')
    fo.write('d,ut,uvw=dget (\''+vis+'\','+str(tel1)+','+str(tel2)+')\n')
    fo.write('np.save(\'d\',d)\n')
    fo.write('np.save(\'ut\',ut)\n')
    fo.write('np.save(\'uvw\',uvw)\n')
    fo.close()
    os.system(CASAPY+' --nologger -c temp.py')
    d,ut,uvw = np.load('d.npy'),np.load('ut.npy'),np.load('uvw.npy')
    os.system('rm d.npy; rm ut.npy; rm uvw.npy')
    return d,ut,uvw

def norm(a,isred=True):
    a = a%(2*np.pi) if isred else a
    np.putmask(a,a>np.pi,a-2.*np.pi)
    np.putmask(a,a<-np.pi,a+2.*np.pi)
    return a

def getap (d,pol=0):
    ph = np.sum(d[pol],axis=0)
    return np.sum(abs(d[pol]),axis=0)/d.shape[1],np.arctan2(ph.imag,ph.real)


def get_uvw_table (t):
    for u in t.select('UVW'):
        try:
            uvw = np.vstack ((uvw,u['UVW']))
        except:
            uvw = np.copy(u['UVW'])
    return uvw

def data_extract (vis):
    global uvw01,uvw02,uvw12,cp012,a01,a02,a12,ra,dec
    chw = taql_calc(vis,'SPECTRAL_WINDOW','CHAN_WIDTH','mean')
    ch0 = taql_calc(vis,'SPECTRAL_WINDOW','REF_FREQUENCY','mean')
    nchan = taql_calc(vis,'SPECTRAL_WINDOW','NUM_CHAN','mean')
    sra,sdec = taql_from(vis,'FIELD','PHASE_DIR')
    nspw = taql_num(vis,'SPECTRAL_WINDOW','NUM_CHAN')
    sra = np.asarray(sra.replace('h',' ').replace('m',' ').split(),dtype='f')
    sdec = np.asarray(sdec.replace('d',' ').replace('m',' ').split(),dtype='f')
    ra = 15.0*(sra[0]+sra[1]/60.0+sra[2]/3600.0)
    dec = sdec[0]+sdec[1]/60.0+sdec[2]/3600.0
    wlength = np.mean(LIGHT/(ch0 + chw*np.arange(nspw*nchan)))
    itel = np.array(idx_tels.get_idx_tels (vis,trname))
    # order of telescopes on baselines 0-1, 0-2, 1-2; -1 if in "wrong" order
    otel = 1.-2.*np.asarray([itel[0]>itel[1],itel[0]>itel[2],itel[1]>itel[2]],dtype=float)
    btel = [min(itel[0],itel[1]),max(itel[0],itel[1]),min(itel[0],itel[2]),\
            max(itel[0],itel[2]),min(itel[1],itel[2]),max(itel[1],itel[2])]
    if have_tables:
        d01,ut01,uvw01 = dget_t (vis, btel[0],btel[1])
        d02,ut02,uvw02 = dget_t (vis, btel[2],btel[3])
        d12,ut12,uvw12 = dget_t (vis, btel[4],btel[5])
    else:
        d01,ut01,uvw01 = dget_c (vis, btel[0],btel[1])
        d02,ut02,uvw02 = dget_c (vis, btel[2],btel[3])
        d12,ut12,uvw12 = dget_c (vis, btel[4],btel[5])
    a01,p01 = getap(d01)
    a02,p02 = getap(d02)
    a12,p12 = getap(d12)
    cp012 = norm(otel[0]*p01-otel[1]*p02+otel[2]*p12,False)
    uvw01 /= wlength
    uvw02 /= wlength
    uvw12 /= wlength
    print trname,'-> antenna numbers:',itel
    print 'Baseline lengths: %s-%s: %dkm %s-%s: %dkm %s-%s: %dkm' % \
       (trname[0],trname[1],int(np.sqrt((uvw01[0]**2).sum())*wlength/1000),\
        trname[0],trname[2],int(np.sqrt((uvw02[0]**2).sum())*wlength/1000),\
        trname[1],trname[2],int(np.sqrt((uvw12[0]**2).sum())*wlength/1000))
    os.system('rm -fr cl_tmp*.ms')
    return itel,np.mean(wlength)

def model_extract (model,itel):
    otel = 1.-2.*np.asarray([itel[0]>itel[1],itel[0]>itel[2],itel[1]>itel[2]],dtype=float)
    re01,im01 = uvw2reim (uvw01,model)
    re12,im12 = uvw2reim (uvw12,model)
    re02,im02 = uvw2reim (uvw02,model)
    ph01 = norm(np.arctan2(im01,re01))
    ph02 = norm(np.arctan2(im02,re02))
    ph12 = norm(np.arctan2(im12,re12))
    clph = norm(ph01*otel[0] - ph02*otel[1] + ph12*otel[2])
    return np.hypot(re01,im01),np.hypot(re02,im02),clph,ph01,ph02,ph12


# ==============================================================

def printmodel (model):
    print '\nCpt     X        Y        flux     width     ellip     PA'
    for i in range(len(model)):
        print '%d  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f'%\
              (i,model[i,0],model[i,1],model[i,2],model[i,3],model[i,4],model[i,5])

def write_skymodel (model,outname='iengine.skymodel'):
    global ra,dec
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


def plotmodel (model,ax):
    for i in range(len(model)):
        e = Ellipse((model[i,MX],model[i,MY]), width=model[i,MW],\
              height=model[i,MW]*(1.-model[i,MR]),angle=model[i,MP])
        ax.add_artist(e)
        e.set_clip_box (ax.bbox)
        e.set_alpha(1.0)
        e.set_facecolor('blue')
    plt.draw()
    printmodel (model)
    write_skymodel (model)

def getclosest (x, y, model):
    dist = 1.0E9
    iret = 0
    for i in range(len(model)):
        thisdist = np.hypot(x-model[i,MX],y-model[i,MY])
        if thisdist<dist:
            dist = thisdist
            iret = i
    return iret

def onkeyclick (event):
    global uvw01,uvw02,uvw12,cp012,a01,a02,a12,asiz,model,trname,itel,ax
    ax = plt.subplot(121)
    plt.cla()
    imod = -1
    try:
        midx = int(event.key)
        if midx < len(model):
             model[midx,MX],model[midx,MY] = event.xdata,event.ydata
        if midx == len(model):
             model = np.vstack((model,\
                   np.array([event.xdata,event.ydata,0.5,0.4,0.0,0.0])))
    except:
        pass
    if event.key == 'm':
        printmodel (model)
        return
    if event.key in 'fFwWrRpP':
        idx = getclosest (event.xdata,event.ydata,model)
        if event.key == 'f':
            model[idx,MF] /= 1.1
        if event.key == 'F':
            model[idx,MF] *= 1.1
        if event.key == 'w':
            model[idx,MW] *= 1.1
        if event.key == 'W':
            model[idx,MW] *= 1.1
        if event.key == 'r':
            model[idx,MR] = min(0.00,model[idx,MR]-0.1)
        if event.key == 'R':
            model[idx,MR] = max(0.99,model[idx,MR]+0.1)
        if event.key == 'p':
            model[idx,MP] = min(0.00,model[idx,MP]-0.1)
        if event.key == 'P':
            model[idx,MP] = max(180.,model[idx,MP]+0.1)
        if event.key == 'a':
            asiz /= 1.1
        if event.key == 'A':
            asiz *= 1.1
        if event.key == 'q':
            plt.close()
            return
    plotmodel (model,ax)
    ax.set_xlim(-asiz,asiz)
    ax.set_ylim(-asiz,asiz)
    plt.grid(b=True)
    A01,A02,CLPH,p01,p02,p12 = model_extract (model,itel)
    ax1 = plt.subplot(322)
    plt.cla();plt.plot(a01);plt.plot(A01);plt.grid(b=True)
    ax2 = plt.subplot(324)
    plt.cla();plt.plot(a02);plt.plot(A02);plt.grid(b=True)
    ax3 = plt.subplot(326)
    plt.cla();plt.plot(cp012);plt.plot(CLPH);plt.grid(b=True)
    plt.draw()

def test(vis):
    global model,itel,ax
    itel,wlength = data_extract(vis)
    A01,A02,CLPH,p01,p02,p12 = model_extract (model,itel)
    slide = [0.1,0.01,0.65,0.02]; slide_init=0.94
    fig = plt.figure(1,figsize=PLOTSIZE)
    ax = plt.subplot(121)
    ax.set_xlim(-asiz,asiz); ax.set_ylim(-asiz,asiz)
    plotmodel (model,ax)
    plt.subplot(322)
    plt.plot(a01);plt.plot(A01);plt.grid(b=True)
    plt.subplot(324)
    plt.plot(a02);plt.plot(A02);plt.grid(b=True)
    plt.subplot(326)
    plt.plot(cp012);plt.plot(CLPH);plt.grid(b=True)
    plt.grid()
    plt.draw()
    cid = fig.canvas.mpl_connect('key_press_event',onkeyclick)
    plt.show()

asiz=10.0
test('SIM5.ms')
