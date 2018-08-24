import matplotlib.pyplot as plt 
import numpy as np
import pynbody

snapshot = np.loadtxt('/z/omarioni/snapshots.txt', dtype='string') #SNAPSHOTS
isnap = snapshot[::-1]

afile = np.loadtxt('/z/omarioni/pyprogram/newdata/M31_masscenter_time.dat')
time = afile[:,0]
xcm =  afile[:,1]
ycm =  afile[:,2]
zcm =  afile[:,3]

xfile = np.loadtxt('/z/omarioni/pyprogram/newdata/part_barra.dat')
ID    = xfile[:,0]
tform = xfile[:,1]
rn_z0 = xfile[:,2]


path = '/srv/cosmdatc/clues/B64_WM3_186592/LG/GAS_SFR/4096_Gasoline/'

kk = range(0,len(ID), 100)
for i in range(0,len(ID)):
    if i in kk:
        print i
    for j in range(0,len(time)-1):
        
        if ((time[j] > tform[i]) & (time[j+1] < tform[i])):

            s=pynbody.load(path + str('%s'%isnap[j]) + '/WMAP3.CLUES.HR.00' + str('%s'%isnap[j]))
            
            IDs = s.star['iord']
            
            xstr= (s.star['pos'].in_units('kpc'))[:,0] - xcm[j] 
            ystr= (s.star['pos'].in_units('kpc'))[:,1] - ycm[j]
            zstr= (s.star['pos'].in_units('kpc'))[:,2] - zcm[j]
            rstr = np.sqrt(xstr**2+ystr**2+zstr**2)    
                                        
            mask3, = np.where(IDs == ID[i])
            
#            print i, ID[i]
            
            archivo1 =  open('/z/omarioni/pyprogram/newdata/tform_barra.dat','a')
            archivo1.write(str('%15d'% ID[i]) +'\t'+
                           str('%12.6f'% tform[i]) +'\t'+
                           str('%12.6f'% rn_z0[i]) +'\t'+
                           str('%12.6f'% rstr[mask3]) +'\n')
            archivo1.close()

            

