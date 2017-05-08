import numpy as np
import scipy as sp
import os,time
import matplotlib.pyplot as plt
from scipy import fftpack as ff
import sys

def mkgauss (naxes,pos,flux,fwhm,axrat=1.0,angle=0.0,ignore=4.0,dodist=False):
# note that total flux = peak flux in a pixel * 1.1331*FWHM**2
# angle is major axis East of North
    a = np.zeros (naxes[0]*naxes[1]).reshape(naxes[1],naxes[0])
    fwhm /= 1.66667
    if axrat==1.0 and angle==0.0:
        for i in range (naxes[1]):
            ydist=float(i)-pos[1]
            for j in range (naxes[0]):
                xdist=float(j)-pos[0]
                if xdist*xdist+ydist*ydist>ignore*ignore*fwhm*fwhm:
                    continue
                if not dodist:
                    a[i,j] = flux*np.exp(-(xdist*xdist+ydist*ydist)/ \
                                    (fwhm*fwhm))/(fwhm*fwhm*np.pi)
                else:
                    a[i,j] = np.hypot(xdist,ydist)
        return a
    sinth = np.sin(angle*np.pi/180.0)
    costh = np.cos(angle*np.pi/180.0)
    r = np.array([-sinth,costh,-costh,-sinth])
    rt = np.array([-sinth,-costh,costh,-sinth])
    sig = np.array([fwhm,0.0,0.0,fwhm*axrat])
    scr1 = mxmul (sig,r)
    scr2 = mxmul (rt, scr1)
    scr1 = mxinv (scr2)
    for i in range(naxes[1]):
        ydist=float(i)-pos[1]
        if abs(ydist)>ignore*fwhm:
            continue
        for j in range (naxes[0]):
            xdist = float(j) - pos[0]
            if abs(xdist)>ignore*fwhm:
                continue
            ex = scr1[0]*xdist+scr1[1]*ydist
            ey = scr1[2]*xdist+scr1[3]*ydist
            if not dodist:
                a[i,j] = (flux/axrat)*np.exp(-(ex*ex+ey*ey))/(fwhm*fwhm*np.pi)
            else:
                a[i,j] = np.hypot(ex,ey)/1.6666667

    return a

def mxmul(a,b):
    output=np.zeros(4)
    output[0]=a[0]*b[0]+a[1]*b[2]
    output[1]=a[0]*b[1]+a[1]*b[3]
    output[2]=a[2]*b[0]+a[3]*b[2]
    output[3]=a[2]*b[1]+a[3]*b[3]
    return output

def mxinv(a):
    det=a[0]*a[3]-a[1]*a[2]
    output=np.array([a[3],-a[1],-a[2],a[0]])/det
    return output

