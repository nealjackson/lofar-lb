from math import *
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import sys, numpy as np
AIPS.userno = int(sys.argv[1])

def time2timera (ttime):
    tday = int(ttime) ; ttime = 24.*(ttime - float(tday))
    thr = int(ttime) ; ttime = 60.*(ttime - float(thr))
    tmin = int(ttime); ttime = 60.*(ttime - float(tmin))
    tsec = int (ttime)
    tlhr,tlmin,tlsec = thr,tmin,tsec-1
    tuhr,tumin,tusec = thr,tmin,tsec+1
    if tusec>59:
        tumin+=1
        tusec-=60
    if tumin>59:
        tuhr+=1
        tumin-=60
    if tlsec<0:
        tlmin-=1
        tlsec+=60
    if tlmin<0:
        tlhr-=1
        tlmin+=60
    return np.array([tday,tlhr,tlmin,tlsec,tday,tuhr,tumin,tusec])

# syntax: checkregular userno inna incl indisk inseq
inna = sys.argv[2]
incl = sys.argv[3]
indisk = 1 if len(sys.argv)<5 else int(sys.argv[4])
inseq = 1 if len(sys.argv)<6 else int(sys.argv[5])
w=WizAIPSUVData(inna,incl,indisk,inseq)
times = []
nant = float(len(w.antennas))
nbas = int( (nant-1.)*0.5*nant)
for wi in w:
    times.append(wi.time)

# Get array of times from the data
    
utime = np.unique(np.asarray(times,dtype='f'))
for u in utime:
    if times.count(u)!=nbas:
        print u,times.count(u)
        new = time2timera(u)
        try:
            flg = np.vstack((flg,new))
        except:
            flg = np.copy(new)

# Flag all times during which one or more baselines is not present

uvflg = AIPSTask('uvflg')
uvflg.indata = AIPSUVData(inna,incl,indisk,inseq)
for t in flg:
    for i in range(8):
        uvflg.timerang[i+1] = int(t[i])
    uvflg.inp()
    uvflg.go()

# Copy the flagged file, which should now be regular as the irregular
# visibilities are dropped during copying (assumes only one FG table)

uvcop = AIPSTask('uvcop')
uvcop.indata = AIPSUVData(inna,incl,indisk,inseq)
uvcop.outdata = AIPSUVData(inna,'REGUL',indisk,1)
uvcop.flagver = 1
uvcop.go()

