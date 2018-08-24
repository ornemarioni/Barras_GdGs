# Es el mismo que el programa que David pero hecho a mi gusto
import numpy as np

def rbin1(x, nbin):
    
    x_sort = np.sort(x)
    
    n = len(x)
    
    delta = n / nbin
    
    med = np.ndarray(nbin)
    
    for i in range(0,nbin):
        med[i] = np.median(x_sort[i*delta:(i+1)*delta])
        
    return med