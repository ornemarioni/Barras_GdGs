import numpy as np
import pynbody

#------------------------------------------------------------------------
vector2 = ('M31', 'MW', 'M33')

path2 = '/z/omarioni/Barras_GdGs/erebos/pyprogram/_data/'

path = '/srv/cosmdatc/clues/B64_WM3_186592/LG/GAS_SFR/4096_Gasoline/'


for l in range(0,3):
    
    file = np.loadtxt(path2 + str('%s'%vector2[l]) +'_tform_z0.dat')
    # file = np.loadtxt(path2 + 'M33_tform_z0_2.dat')
    ID    = file[:,0]
    tform = file[:,2]

    sort = np.argsort(tform)

    sort_tform = tform[sort][::-1]
    sort_ID    = ID[sort][::-1]

#     file2 = np.loadtxt(path2 + 'M33_masscenter.dat')
    file2 = np.loadtxt(path2 + str('%s'%vector2[l]) +'_masscenter.dat')
    
    time = file2[:,0]
    xcm  = file2[:,1]
    ycm  = file2[:,2]
    zcm  = file2[:,3]
    
    time_aux = np.zeros(len(ID))

    k = 0
    
    for i in range(0,len(ID)):
        for j in range(k, len(time)):
            if sort_tform[i] > time[0]:
                time_aux[i] = time[0]
                break
                
            if sort_tform[i] < time[j]:
                time_aux[i] = time[j]
    #             print time_aux[i]
                if time_aux[i] < time_aux[i-1]:
                    k = k + 1
                break
                
    
    snapshot = np.loadtxt('/z/omarioni/snapshots.txt', dtype='string') #SNAPSHOTS
    isnap = snapshot[::-1]

#     archivo = open(path2 + 'M33_tform_particles.dat', 'a')
    archivo = open(path2 + str('%s'%vector2[l]) +'_tform_particles.dat', 'a')
    
    time_aux = time_aux[::-1]
    
    
    for i in range(0, len(time)):

        mask, = np.where(time_aux == time[i])

        if len(mask) == 0:
            continue

        s = pynbody.load(path + str('%s'%isnap[i])+'/WMAP3.CLUES.HR.00'+ str('%s'%isnap[i]))

        pstr = s.star['pos'].in_units('kpc')
        IDs  = s.star['iord']
        tf   = s.star['tform'].in_units('Gyr')

        particles = np.isin(IDs, np.int_(sort_ID[mask]))

        xstr = (pstr[:,0][particles]-xcm[i])
        ystr = (pstr[:,1][particles]-ycm[i])
        zstr = (pstr[:,2][particles]-zcm[i])

        rstr = np.sqrt(xstr**2 + ystr**2 + zstr**2)

        data = np.ndarray([len(mask), 3])
        data[:,0] = IDs[particles]
        data[:,1] = rstr
        data[:,2] = tf[particles]

        np.savetxt(archivo, data, fmt=('%12d','%12.6f','%12.6f'))

    archivo.close()