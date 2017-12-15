import numpy as np
import bines2 as bines 

def a2max(m,x,y,nbin):
    
    Rcil = np.sqrt(x**2 + y**2)
    rbin = bines.rbin1(Rcil, nbin)
    
    A2v  = np.ndarray(nbin)
    phiv = np.ndarray(nbin)
    
    for j in range(0,nbin-1):
        mask, = np.where((Rcil > rbin[j]) & (Rcil <= rbin[j+1]))
        
        xn = x[mask]
        yn = y[mask]
        
        titaj = np.arctan2(yn,xn)
        
        a0 = sum(m[mask])
        a2 = sum(m[mask]*np.cos(2.*titaj))
        b2 = sum(m[mask]*np.sin(2.*titaj))
        
        A2v[j]  = np.sqrt(a2**2 + b2**2) / a0
        phiv[j] = np.arctan2(b2,a2) / 2.
        
        
    A2max_arg = np.argmax(A2v)
    
    A2max   = max(A2v)
    rbinmax = rbin[A2max_arg]
    phimax  = phiv[A2max_arg]
    
    return A2max, rbinmax, phimax



def a2(m,x,y,nbin):
    
    Rcil = np.sqrt(x**2 + y**2)
    rbin = bines.rbin1(Rcil, nbin)
    
    A2v  = np.ndarray(nbin)
    phiv = np.ndarray(nbin)
    
    for j in range(0,nbin-1):
        mask, = np.where((Rcil > rbin[j]) & (Rcil <= rbin[j+1]))
        
        xn = x[mask]
        yn = y[mask]
        
        titaj = np.arctan2(yn,xn)
        
        a0 = sum(m[mask])
        a2 = sum(m[mask]*np.cos(2.*titaj))
        b2 = sum(m[mask]*np.sin(2.*titaj))
        
        A2v[j]  = np.sqrt(a2**2 + b2**2) / a0
        phiv[j] = np.arctan2(b2,a2) / 2.
        
    
    return A2v, phiv, rbin




