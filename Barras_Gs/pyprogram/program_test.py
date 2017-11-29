import numpy as np
import matplotlib.pyplot as plt
import pynbody
snap = np.loadtxt('/z/dalgorry/snapshots.txt', dtype='string')
snap = sorted(snap, reverse=True)

halo1 = np.loadtxt('/z/dalgorry/datos_salida2/merger_tree-h1.dat')
v = halo1[:,2]
v = np.int_(v)

fo=open('/z/dalgorry/datos_salida2/time_snap_test.dat','w')

for j in range(0,len(snap)-1,1):
    
    s=pynbody.load('/srv/cosmdatc/clues/B64_WM3_186592/LG/GAS_SFR/4096_Gasoline/'+str('%s'%snap[j])+'/WMAP3.CLUES.HR.00'+str('%s'%snap[j]))
    h = s.halos()
    h1 = h[v[j]]
    
    time = s.properties['time'].in_units('Gyr')
    aexp = s.properties['a']
    rvir = h1.properties['Rvir']
    hh = s.properties['h']
    
    cen_pot = pynbody.analysis.halo.center(h1,mode='pot',retcen=True).in_units('kpc')
    x = h1.dm['pos'][:,0].in_units('kpc')-cen_pot[0]
    y = h1.dm['pos'][:,1].in_units('kpc')-cen_pot[1]
    z = h1.dm['pos'][:,2].in_units('kpc')-cen_pot[2]
    r = np.sqrt(x**2+y**2+z**2)
    
    r_sort = np.sort(r)
    
    fo.write(str('%s'% snap[j])+'\t'+
            str('%12.6f'% time)+'\t'+
            str('%12.6f'% rvir)+'\t'+
            str('%12.6f'% r_sort[-1])+'\t'+
            str('%12.6f'% aexp)+'\t'+
            str('%12.6f'% hh)+'\n')
    fo.flush()
    print j
