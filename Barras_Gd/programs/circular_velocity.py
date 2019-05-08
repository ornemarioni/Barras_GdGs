import numpy as np
import bines2

#Unidades del G [kpc Msol^-1 (km/s)^2]

def Vcirc(r, m, rmax=200.,G = 4.299e-6):

    
    limit = np.where(r < rmax)
    r_ind = np.argsort(r[limit])

    M = np.cumsum(m[limit][r_ind])

    Vcir = np.sqrt(G * M / r[limit][r_ind])

    return Vcir, r[limit][r_ind]


#con bineado
def Vc_bin(r, m, nbin,G = 4.299e-6):

    med, nodos = bines2.rbin1(r, nbin)
    
    Vcir = np.zeros(nbin)
    
    for i in range(nbin):
        limit = np.where(r < med[i])

        M = np.sum(m[limit])

        Vcir[i] = np.sqrt(G * M / med[i])

    return Vcir, med