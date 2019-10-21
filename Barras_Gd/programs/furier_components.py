import numpy as np
import bines2 as bines 

def am(m,x,y,nbin, modo):
    rcil = np.sqrt(x**2 + y**2)
    rbin, nodos = bines.rbin1(rcil, nbin)
    
    Amv  = np.ndarray(nbin)
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
        am = np.sum(m[mask]*np.cos(modo*titaj))
        bm = np.sum(m[mask]*np.sin(modo*titaj))
        
        Amv[j]  = np.sqrt(am**2 + bm**2) / a0
        phiv[j] = np.arctan2(bm,am) / 2.
        

    return A2v, phiv, rbin


