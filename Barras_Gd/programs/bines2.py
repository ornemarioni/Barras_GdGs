# Es el mismo que el programa que David pero hecho a mi gusto
import numpy as np

def rbin1(x, nbin):
    
    x_sort = np.sort(x)
    
    n = len(x)
    
    delta = n / nbin
    
    med = np.zeros(nbin)
    
    nodos    = np.zeros(nbin+1)
    nodos[0] = x_sort[0]
    
    for i in range(0,nbin):
        med[i]     = np.median(x_sort[i*delta:(i+1)*delta])
        nodos[i+1] = x_sort[i*delta:(i+1)*delta][-1]
    
    return med, nodos


###aca se agregaron los nodos,tengo que modificar los programas que usan esto 