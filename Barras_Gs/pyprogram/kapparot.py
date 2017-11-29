import numpy as np

def kapa(m,x,y,z,vx,vy,vz):
        rcil2=(x**2+y**2)
        kk=np.where(rcil2 !=0)

        ekin=sum(m[kk]*(vx[kk]*vx[kk]+vy[kk]*vy[kk]+vz[kk]*vz[kk]))  

        jz=x[kk]*vy[kk]-y[kk]*vx[kk]
        jz1=m[kk]*jz**2/rcil2[kk]
        jzacum=sum(jz1)

        kapparot=jzacum/ekin

        return kapparot
         

        
