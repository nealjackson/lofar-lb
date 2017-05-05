import numpy as np,matplotlib; from matplotlib import pyplot as plt
import scipy; from scipy import special

def pattern7 (cen, throw)
    p6 = np.array(np.arange(0,2.*np.pi,np.pi/3))
    return np.vstack((np.array([cen[0]+throw*np.cos(p6),cen[1]+throw*np.sin(p6)]).T,cen))

def pattern19 (cen, throw):
    p7a = pattern7 (cen,throw)[:-1]
    p7b = pattern7 (cen, 2.*throw)[:-1]
    p7c = pattern7 (cen, np.sqrt(3.)*throw)[:-1] - cen
    xnew = 0.5*np.sqrt(3)*p7c[:,0]+0.5*p7c[:,1] + cen[0]
    ynew = -0.5*p7c[:,0]+0.5*np.sqrt(3)*p7c[:,1] + cen[1]
    p7c = np.dstack((xnew,ynew))[0]
    return np.vstack((p7a,p7b,p7c,cen))

def mkoffsets (ra, dec):   # ra,dec in degrees
    a = pattern19(np.array([0,0]),1600./3600.)    # these values optimized
    b = pattern7(np.array([0,0]),600./3600.)
    apos,bpos = np.zeros((19,2)),np.zeros((19,7,2))
    apos[:,0] = ra+a[:,0]/np.cos(np.deg2rad(dec))
    apos[:,1] = dec+a[:,1]
    for i in range(19):
        bpos[i,:,0] = apos[i,0] + b[:,0]/np.cos(np.deg2rad(apos[i,1]))
        bpos[i,:,1] = apos[i,1] + b[:,1]
    return apos,bpos
