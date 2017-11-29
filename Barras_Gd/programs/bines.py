import numpy as np

def rbin1(x,y,nbin):
    k  = np.argsort(x)
    xx = x[k]
    yy = y[k]
    n  = np.shape(xx)[0] 
    delta = n/nbin
    med = np.ndarray(nbin)

    for i in np.arange(1, nbin+1, 1):
           med[i-1] = np.median(xx[((i-1)*delta):(i*delta)])
 
    return med
