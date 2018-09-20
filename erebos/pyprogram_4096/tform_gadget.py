import numpy as np
import gadget_reader as gd
import time_conversion as tiempo
from joblib import Parallel, delayed

# Parametros cosmologicos
#--------------------------------------------------------------------------
h = 0.732
G = 4.299e-6
a0=1.
H0 = h*100
omega_lambda=0.716
omega_matter=0.1277/(h**2.)
omega0 = omega_lambda + omega_matter
#------------------------------------------------------------------------
vector2 = ('M31', 'MW', 'M33')

path2 = '/z/omarioni/Barras_GdGs/erebos/pyprogram_4096/'

aa = np.loadtxt(path2 + 'redshift_outputs.txt')
aexp  = aa[:,2]

path  = '/srv/cosmdatc/clues/B64_WM3_186592/LG/GAS_SFR/4096/SNAPS/snap_'


def sarasa(l):
# for l in range(0,3):
    
    file = np.loadtxt(path2 + str('%s'%vector2[l]) +'_tform_z0.dat')
    ID    = file[:,0]
    tform = file[:,2]
    
    sort = np.argsort(tform)
    
    sort_tform = tform[sort][::-1]
    sort_ID    = ID[sort][::-1]
    
    file2 = np.loadtxt(path2 +  str('%s'%vector2[l]) +'_masscenter.dat')
    time = file2[:,0]
    xcm  = file2[:,1]
    ycm  = file2[:,2]
    zcm  = file2[:,3]
    
    
    time_aux = np.zeros(len(ID))

    k = 1
    
    for i in range(0,len(ID)):
        if sort_tform[i] > time[0]:
            time_aux[i] = time[0]
            continue
            
        if sort_tform[i] < time[-1]:
            time_aux[i] = time[-1]
            continue
            
        for j in range(k, len(time)):
            
            if  sort_tform[i] > time[j]:
                time_aux[i] = time[j-1]
    #             print time_aux[i]
                if time_aux[i] < time_aux[i-1]:
                    k = k + 1
                break

    isnap = np.arange(496,0,-1)
    
    archivo = open(path2 + str('%s'%vector2[l]) + '_tform_particles.dat', 'a')

    for i in range(0, len(time)):

        mask, = np.where(time_aux == time[i])

        if len(mask) == 0:
            break

        sim = gd.Open(path + str('%03d'%isnap[i]), endian='Big', gadget_type=2, verbose=False)

        pstr = sim.Read('POS ',4)
        IDs  = sim.Read('ID  ',4)
        tf   = sim.Read('AGE ',4)

        z   = a0/tf - 1.
        Ht  = H0*np.sqrt(omega_lambda+(1-omega0)*(1+z)**2+omega_matter*(1+z)**3)
        tf2 = tiempo.conv(z, h, omega_lambda, omega_matter)


        particles = np.isin(IDs, sort_ID[mask])


        xstr = (pstr[0,:][particles]-xcm[i])*aexp[isnap[i]]/h
        ystr = (pstr[1,:][particles]-ycm[i])*aexp[isnap[i]]/h
        zstr = (pstr[2,:][particles]-zcm[i])*aexp[isnap[i]]/h

        rstr = np.sqrt(xstr**2 + ystr**2 + zstr**2)

        data = np.ndarray([len(mask), 3])
        data[:,0] = IDs[particles]
        data[:,1] = rstr
        data[:,2] = tf2[particles]
        
        np.savetxt(archivo, data, fmt=('%15d','%12.6f','%12.6f'))
        
    archivo.close()
    

with Parallel(n_jobs=3, prefer="threads") as par:
    par(delayed(sarasa)(ll)for ll in range(0,3))

  