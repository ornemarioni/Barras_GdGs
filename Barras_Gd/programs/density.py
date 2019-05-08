import numpy as np
import bines2

#----------------------------------------------------------------------
# Densidad superficial de masa
#----------------------------------------------------------------------
def surf_density(R,m,nbin):
    
    med, nodos = bines2.rbin1(R, nbin)
    
    sigma = np.zeros(nbin)
    
    for i in range(nbin):
        area = np.pi*(nodos[i+1]**2 - nodos[i]**2)
        
        mask, = np.where((R < nodos[i+1]) & (R > nodos[i]))
        
        sigma[i] = np.sum(m[mask])/area
        
    return sigma, med
        
    
#----------------------------------------------------------------------
# Densidad volumetrica de masa
#----------------------------------------------------------------------    
def vol_density(r, m, nbin):
    
    med, nodos = bines2.rbin1(r, nbin)
    
    rho = np.zeros(nbin)
    
    for i in range(nbin):
        vol = (4/3.)*np.pi*(nodos[i+1]**3 - nodos[i]**3)
        
        mask, = np.where((r < nodos[i+1]) & (r > nodos[i]))
        
        rho[i] = np.sum(m[mask])/vol
        
    return rho, med