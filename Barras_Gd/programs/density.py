import numpy as np
import bines2

#----------------------------------------------------------------------
# Densidad superficial de masa con bines equal number
#----------------------------------------------------------------------
def surf_density(R,m,nbin):
    
    if len(R) != len(m):
        print 'Error: vector lengths do not match'
        
    med, nodos = bines2.rbin1(R, nbin)
    
    sigma = np.zeros(nbin)
    
    for i in range(nbin):
        area = np.pi*np.abs(nodos[i+1]**2 - nodos[i]**2)
        
        mask, = np.where((R < nodos[i+1]) & (R > nodos[i]))
        
        sigma[i] = np.sum(m[mask])/area
        
    return sigma, med

#----------------------------------------------------------------------
# Densidad superficial de masa con bines equiespaciados
#----------------------------------------------------------------------
def surf_density2(R,m,nbin):
    
    if len(R) != len(m):
        print 'Error: vector lengths do not match'
        
    nodos = np.linspace(R.min(),R.max(),nbin+1)
    
    width = (R.max()-R.min())/nbin
    
    med = np.linspace(R.min()+width/2., R.max()-width/2., nbin)
    
    sigma = np.zeros(nbin)
    
    for i in range(nbin):
        area = np.pi*np.abs(nodos[i+1]**2 - nodos[i]**2)
        
        mask, = np.where((R < nodos[i+1]) & (R > nodos[i]))
        
        sigma[i] = np.sum(m[mask])/area
        
    return sigma, med

#----------------------------------------------------------------------
# Densidad superficial de masa con bines equiespaciados (en x o y)
#----------------------------------------------------------------------
def surf_density_x(x,y,m,nbin):
    
    if len(x) != len(m):
        print 'Error: vector lengths do not match'
        
    nodos_ = np.linspace(0,x.max(),nbin+1)
    nodos  = np.concatenate([-nodos_[::-1],nodos_])
    
    width = x.max()/nbin
    
    med_ = np.linspace(width/2., x.max()-width/2., nbin)
    med  = np.concatenate([-med_[::-1],med_])
    
    sigma = np.zeros(nbin*2)
    
    for i in range(nbin*2):
        area = width**2
        
        mask, = np.where((x < nodos[i+1]) & (x > nodos[i]) & (y<width/2.) & (y>-width/2.))
        
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