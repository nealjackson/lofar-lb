import numpy as np, os, sys, glob, astropy
from lofipi_aips import *
import Wizardry; from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io.fits import getdata,getheader
from correlate import *
indisk = 1
inseq = 1
inlog='temp.log'
AIPS.userno = 340
BAND = 'C'
pos_3C48 = [24.42208,33.15976]
pos_3C286 = [202.7845,30.50915]
pos_3C138 = [80.2911917,16.63946]
pos_3C147 = [85.650574,49.852008]
FITSDIR, TABLESDIR = '../newdata/', './newresults/'
CASA = '/mirror3/scratch/scratch/casa-release-5.4.1-32.el7/bin/casa'
IMSIZ0 = 1024
IMSF = 128
CELL = 0.1
NITER = 800
ILIM=12

def logprint (inna, text):
    fo = open(inlog,'a')
    fo.write ('proc> '+inna+': '+text+'\n')
    fo.flush()
    fo.close()
    print (text)

def exists (inna, incl, indi, inseq):
    for i in AIPSCat()[indi]:
        if i['name']==inna and i['klass']==incl and i['seq']==inseq:
            return True
    return False

def zapit (inna,incl,indi,inseq):
    for i in AIPSCat()[indi]:
        if inseq==0:
            if i['name']==inna and i['klass']==incl:
                if i['type']=='UV':
                    logprint (inna, 'Deleting %s.%s.%d %d'%(inna,incl,inseq,indi))
                    AIPSUVData(inna,incl,indi,inseq).clrstat()
                    AIPSUVData(inna,incl,indi,inseq).zap()
                else:
                    logprint (inna, 'Deleting %s.%s.%d %d'%(inna,incl,inseq,indi))
                    AIPSImage(inna,incl,indi,inseq).clrstat()
                    AIPSImage(inna,incl,indi,inseq).zap()
            break   #  if inseq=0 only do it once
        else:
            if i['name']==inna and i['klass']==incl and i['seq']==inseq:
                if i['type']=='UV':
                    AIPSUVData(inna,incl,indi,inseq).clrstat()
                    AIPSUVData(inna,incl,indi,inseq).zap()
                else:
                    AIPSImage(inna,incl,indi,inseq).clrstat()
                    AIPSImage(inna,incl,indi,inseq).zap()

def globalzapit (indisk):
    for i in AIPSCat()[indisk]:
        zapit (i['name'],i['klass'],indisk,i['seq'])

def getsrc (inna, indisk, incl):
    print (inna, incl, indisk)
    data = WizAIPSUVData(inna, incl, indisk, 1)
    nchan = data[0].visibility.shape[1]
    su = data.table('SU',1)
    sources,ra,dec,intent = [],[],[],[]
    have_fluxcal,have_targlist = False,False
    if os.path.isfile('targets_list'):
        have_targlist = True
        targlist = np.loadtxt('targets_list',dtype='str')
    for i in range(len(su)):
        sources.append(su[i].source.rstrip())
        ra.append(su[i].raobs if su[i].raobs>0.0 else 360+su[i].raobs)
        dec.append(su[i].decobs)
        if np.hypot(ra[i]-pos_3C48[0],dec[i]-pos_3C48[1]) < 0.01:
            intent.append('3C48')
            have_fluxcal = True
        elif np.hypot(ra[i]-pos_3C286[0],dec[i]-pos_3C286[1]) < 0.01:
            intent.append('3C286')
            have_fluxcal = True
        elif np.hypot(ra[i]-pos_3C138[0],dec[i]-pos_3C138[1]) < 0.01:
            intent.append('3C138')
            have_fluxcal = True
        elif np.hypot(ra[i]-pos_3C147[0],dec[i]-pos_3C147[1]) < 0.01:
            intent.append('3C147')
            have_fluxcal = True
        else:
            if have_targlist:
                intent.append('T' if sources[i] in targlist else 'C')
            else:
                logprint (inna, 'No targets_list file found, exiting')
                exit(0)
    if not have_fluxcal:
        logprint (inna, 'No flux calibrator found, exiting')
        exit(0)
    return np.asarray(sources,dtype='str'), \
           np.asarray(ra, dtype='float'), \
           np.asarray(dec, dtype='float'), \
           np.asarray(intent, dtype='str'), nchan

def string2casa (coord):
    coord=coord.replace('h',':')
    coord=coord.replace('s','')
    coord=coord.replace('d','.')
    coord=coord.replace('m',':',1)
    coord=coord.replace('m','.')
    return coord

def getfields (ra,dec,sourcename):
    nfield = 1
    try:
        first = np.loadtxt('first_simple.txt')
    except:
        os.system('wget http://www.jb.man.ac.uk/~njj/first_simple.txt')
        first = np.loadtxt('first_simple.txt')
    acor = correlate(np.array([[ra,dec]]),0,1,first,0,1,0.5)
    logprint (sourcename, 'Identified %d nearby sources'%int(len(acor)))
    coord = SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
    phase = ['J2000 %s'%string2casa(coord.to_string('hmsdms'))]
    for aline in acor:
        first_index = int(aline[1])
        thisdist, thisflux = aline[2]*3600, first[first_index,3]
        thisra, thisdec = first[first_index,0], first[first_index,1]
        thiscoord = SkyCoord(thisra*u.deg,thisdec*u.deg,frame='icrs')
        thisstr = thiscoord.to_string('hmsdms')
        halfmap = 0.5*IMSIZ0*CELL/3600.0
#   ignore if the first source is weak or far away or both
        if thisdist<200. or (thisdist<500. and thisflux>5.) or (thisdist<1000. and thisflux>20.):
            pass
        else:
            logprint (sourcename, 'Rejecting source %s too faint/far away'%(thisstr))
            continue
#   ignore if the first source is already in the main field
        if abs(dec-thisdec)<halfmap and \
           abs((ra-thisra)*np.cos(np.deg2rad(dec)))<halfmap:
            logprint (sourcename, 'Rejecting source %s already in main field'%(thisstr))
            continue
        coord = SkyCoord(first[first_index,0]*u.deg,first[first_index,1]*u.deg)
        phase.append('J2000 %s'%string2casa(coord.to_string('hmsdms')))
        nfield += 1
        logprint (sourcename, 'Total of %d sources in field for imaging'%nfield)
    return nfield, phase

def get_unflag_idx (val):
    r = np.asarray(np.arange(21)*len(val)/20,dtype='int')
    refval = np.median(val[r[19]:])
    idx = 0
    for i in range(18,0,-1):
        thisval = np.median(val[r[i]:r[i+1]])
        if thisval < 0.9*refval:
             idx = r[min(i+2,19)]
             break
    return idx

def unflag_fluxcal (inna, incl, indisk, fluxno):
    # find index number of fluxcal in data
    data = WizAIPSUVData (inna, incl, indisk, 1)
    start_idx, end_idx = -1,0
    npos, nall = 0,0
    val = np.array([],dtype='float')
    logprint (inna, 'Starting pass to identify data')
    for i in range(len(data)):
        if data[i].source == fluxno and start_idx==-1:
            start_idx = end_idx = i
        if start_idx >-1 and data[i].source == fluxno:
            this = np.ravel(data[i].visibility[:,:,:,2])
            real = np.ravel(data[i].visibility[:,:,0,0])
            imag = np.ravel(data[i].visibility[:,:,0,1])
            amp = np.hypot(real,imag)
            val = np.append(val,np.median(amp))
            end_idx = i
            npos += len(this[this>0])
            nall += len(this)
        if start_idx >-1 and data[i].source != fluxno:
            break
    logprint (inna, 'Fraction of good data: %d/%d'%(npos,nall))
    # brute-force unflagging by setting the weights
    if float(npos)/float(nall)<0.2:
        unflag_idx = get_unflag_idx (val)
        data = WizAIPSUVData (inna, incl, indisk, 1)
        logprint (inna, 'Unflagging fluxcal from indices %d-%d'%\
                (start_idx+unflag_idx,end_idx))
        for visibility in data[start_idx+unflag_idx:end_idx]:
            visibility.visibility[:,:,:,2]=abs(visibility.visibility[:,:,:,2])
            visibility.update()
        return True
    else:
        return False



def make_casastr (nfield, source, tc):
    if nfield==1:
        imname_str = '\'%s\'' % source
        phase_str = '\'%s\'' % tc[0]
        imsize_str = '[%d,%d]'%(IMSIZ0,IMSIZ0)
    else:
        imname_str = '[\''+source+'\''
        phase_str = '['
        imsize_str = '[[%d,%d]'%(IMSIZ0,IMSIZ0)
        for j in range(nfield-1):
            phase_str+=('\''+tc[j]+'\',')
        phase_str+=('\''+tc[nfield-1]+'\']')
        for j in range(1,nfield-1):
            imname_str=imname_str+',\'%s_%02d\''%(source,j)
            imsize_str+=',[%d,%d]'%(IMSF,IMSF)
        imname_str=imname_str+(',\'%s_%02d\']'%(source,nfield-1))
        imsize_str+=',[%d,%d]]'%(IMSF,IMSF)
    return imname_str,phase_str,imsize_str


def proc_casaimg (fitsfile,tasav_file):
    globalzapit (indisk)
    inna = fitsfile.split('/')[-1].split('F.fits')[0]
    if not exists (inna, 'UVDATA', indisk, inseq):
        fitld (fitsfile, inna, indisk, 'UVDATA', doindxr=1, iplog=inlog)
    
    sources, ra, dec, intent, nchan = getsrc (inna, indisk, 'UVDATA')
    for i in range (len(intent)):
        logprint(inna,'source %s, intent %s'%(source[i],intent[i]))
        if intent[i][:2]=='3C':
            fluxc = str(intent[i])
            
    autoquack = unflag_fluxcal (inna, 'UVDATA', indisk, 1+np.where(intent==fluxc)[0][0])
    
    if not exists (inna, 'TABLES', indisk, inseq):
        fitld (tasav_file, inna, indisk, 'TABLES', doindxr=1, iplog=inlog)
    
    n_fluxc = sources[intent==fluxc]
    n_target = sources[intent=='T']
    n_phcal = np.append(sources[intent=='C'],sources[intent=='fluxc'])
    if not exists (fluxc+'_'+BAND, 'MODEL', indisk, 1):
        calrd (fluxc, BAND, indisk, iplog=inlog)
    
    tacop (inna, indisk, 'TABLES', inna, indisk, 'UVDATA', 'SU', 1, iplog=inlog)
    tacop (inna, indisk, 'TABLES', inna, indisk, 'UVDATA', 'CL', 2, iplog=inlog)
    bpass (inna, 'UVDATA', indisk, docalib=1, calsou=n_fluxc.tolist(), iplog=inlog)
    for i in range(len(sources)):
        if intent[i]!='T':
            continue
        split (inna, 'UVDATA', sources[i], indisk)
        clip (sources[i][:ILIM], 'SPLIT', indisk=1, aclip=0.8)
        nfield, tc = getfields (ra[i], dec[i], sources[i])
        os.system('rm ./'+sources[i]+'_cal.fits')
        if nfield==1:
            avspc (sources[i][:ILIM], 'UVCOP', 1, 'UAVSP', 4)
            uvavg (sources[i][:ILIM], 'UAVSP', 1, 'UAVG', yinc=8, zinc=2)
            fittp (sources[i][:ILIM], 'UAVG', indisk, 1, './'+sources[i]+'_cal.fits')
        else:
            fittp (sources[i][:ILIM], 'UVCOP', indisk, 1, './'+sources[i]+'_cal.fits')
        
        os.system('rm -fr %s_cal.ms'%(sources[i]))
        fo = open('casacommands.txt','w')
        fo.write('import os\n')
        fo.write('importuvfits(fitsfile=\''+sources[i]+'_cal.fits\',vis=\''+sources[i]+'_cal.ms\')\n')
        fo.write('statwt(\''+sources[i]+'_cal.ms\',datacolumn=\'DATA\')\n')
#  only make the natural maps
#        imname_str,phase_str,imsize_str = make_casastr (nfield,sources[i]+'B',tc)
#        fo.write('clean(vis=\'%s_cal.ms\',imsize=%s,cell=\'%.1farcsec\',phasecenter=%s,mode=\'mfs\',niter=%d,weighting=\'briggs\',imagename=%s)\n'%(sources[i],imsize_str,CELL,phase_str,NITER,imname_str))
#        fo.write('os.system(\'rm %s_cbri.fits\')\n'%(sources[i]))
#        fo.write('exportfits(imagename=\'%sB.image\',fitsimage=\'%s_cbri.fits\')\n'%(sources[i],sources[i]))
        imname_str,phase_str,imsize_str = make_casastr (nfield,sources[i]+'N',tc)
        fo.write('clean(vis=\'%s_cal.ms\',imsize=%s,cell=\'%.1farcsec\',phasecenter=%s,mode=\'mfs\',niter=%d,weighting=\'natural\',imagename=%s)\n'%(sources[i],imsize_str,CELL,phase_str,NITER,imname_str))
        fo.write('os.system(\'rm %s_cnat.fits\')\n'%(sources[i]))
        fo.write('exportfits(imagename=\'%sN.image\',fitsimage=\'%s_cnat.fits\')\n'%(sources[i],sources[i]))
        fo.close()
        os.system(CASA+' -c casacommands.txt')
        os.system('mv %s* casaresults'%(sources[i]))
        os.system('rm *.last')
        os.system('mv temp.log casaresults/%s.log'%(sources[i]))
        os.system('casaresults/tidy')
    

# Done so far: 38766441, 38766904, 38768393, 38766614, 38768661, 38771665, 38767180,38768534,38769646,38770337,38771097,38771223,38911390,38911504,38911724'38911947','38912052','38768780','38769410','38770477'

#proc_casaimg (FITSDIR+'%sF.fits'%sys.argv[1],TABLESDIR+'%s_tasav.fits'%sys.argv[1])
#allf = ['38769283','38770606','38770215']
allfiles = np.sort(glob.glob('newresults/4*tasav.fits'))
allf = []
for i in allfiles:
    allf.append(i.split('/')[-1].split('_')[0])
for i in allf[7:]:
    try:
        proc_casaimg(FITSDIR+'%sF.fits'%i,TABLESDIR+'%s_tasav.fits'%i)
    except:
        pass
