#Parametro de Toomre
import numpy as np
import density as den

def toomre(r, m, Vr, nbin, G = 4.299e-6)

    med, nodos = bines2.rbin1(r, nbin)
    
    Vcir = np.zeros(nbin)
    M    = np.zeros(nbin)
    
    for i in range(nbin):
        limit = np.where(r < med[i])
        M[i]  = np.sum(m[limit])
        Vcir[i] = np.sqrt(G * M[i] / med[i])
        

        mask, = np.where((nodos[i] < R) & (nodos[i+1] > R))
        n = len(mask)      
        Vr_mean = np.mean(Vr[mask])
        sigma_R[i] = np.sqrt(sum((Vr[mask] - Vr_mean)**2)/n)

    
    surf_den, rden = den.surf_density(r, m, nbin)
    
    a0 = 2.*Vcir/med**2.
    a1 = med/2.
    a2 = np.diff(M)/np.diff(med)
    a3 = (a2/M**(3./2)) - (1./med)
    
    k2 = a0 * (1. + a1*a3)
    k  = np.sqrt(k2)
    
    Q = sigma_R * k /(3.36 * G * sigma_surf)
    
    return Q, med