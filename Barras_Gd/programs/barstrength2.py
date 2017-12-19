import numpy as np
import bines2 as bines 

def a2max(m,x,y,nbin):
    rcil = np.sqrt(x**2 + y**2)
    rbin = bines.rbin1(rcil, nbin)
    
    A2v  = np.ndarray(nbin)
    phiv = np.ndarray(nbin)
    
    n = len(x)
    delta = n / nbin
    
    r_sort = np.sort(rcil)
    
    for j in range(0, nbin):
        mask, = np.where((rcil > r_sort[j*(delta-1)]) & (rcil <= r_sort[(delta-1)*(j+1)]))
        
        xn = x[mask]
        yn = y[mask]
        
        titaj = np.arctan2(yn,xn)
        
        a0 = np.sum(m[mask])
        a2 = np.sum(m[mask]*np.cos(2.*titaj))
        b2 = np.sum(m[mask]*np.sin(2.*titaj))
        
        A2v[j]  = np.sqrt(a2**2 + b2**2) / a0
        phiv[j] = np.arctan2(b2,a2) / 2.
        
    A2max_arg = np.argmax(A2v)
    
    A2max   = A2v[A2max_arg]
    rbinmax = rbin[A2max_arg]
    phimax  = phiv[A2max_arg]
 
    return A2max, rbinmax, phimax


def a2(m,x,y,nbin):
    rcil = np.sqrt(x**2 + y**2)
    rbin = bines.rbin1(rcil, nbin)
    
    A2v  = np.ndarray(nbin)
    phiv = np.ndarray(nbin)
    
    n = len(x)
    delta = n / nbin
    
    r_sort = np.sort(rcil)
    
    for j in range(0, nbin):
        mask, = np.where((rcil > r_sort[j*(delta-1)]) & (rcil <= r_sort[(delta-1)*(j+1)]))
        
        xn = x[mask]
        yn = y[mask]
        
        titaj = np.arctan2(yn,xn)
        
        a0 = np.sum(m[mask])
        a2 = np.sum(m[mask]*np.cos(2.*titaj))
        b2 = np.sum(m[mask]*np.sin(2.*titaj))
        
        A2v[j]  = np.sqrt(a2**2 + b2**2) / a0
        phiv[j] = np.arctan2(b2,a2) / 2.
        
 
    return A2v, phiv, rbin







