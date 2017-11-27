"""This file calculate the Jcir parameter (Abadi 2003).
e=Energy
jz=Angular momentum
gbin=length of the grid map e vs jz. Default value=200
nbines=length of the wraping curve. Default value=30 
thresholds== limits on density of particles per each bin. Default value=0"""

import numpy as np


def jcirc(e,jz,gbines=200,nbines=30,threshold=0):

	#preparing wrapping curve (e vs jz)

        n=np.size(e)
        emin=min(e)
        emax=0.
        jzmin=min(jz)
        jzmax=max(jz)
        
        
        grid=np.zeros([gbines,gbines])
        
        for i in range(n):	
        	indice1=np.int32((e[i]-emin)/(emax-emin)*(gbines-1))
        	indice2=np.int32((jz[i]-jzmin)/(jzmax-jzmin)*(gbines-1))
        	grid[indice1,indice2]=grid[indice1,indice2]+1			
        
        
        print 'Determining the local particle density...'
        
        rho=np.ndarray([n])
        
        for i in range(n):	
        	indice1=np.int32((e[i]-emin)/(emax-emin)*(gbines-1))
        	indice2=np.int32((jz[i]-jzmin)/(jzmax-jzmin)*(gbines-1))
        	rho[i]=grid[indice1,indice2]			
        
        xx=np.ndarray([nbines])
        yy=np.ndarray([nbines])
        ll=np.array(range(nbines))
        
        xx[:]=emin+(ll[:]+1.0)*(emax-emin)/(nbines-1)
        
        yy[:]=-1.e-10
        emini=-1.e-10
        
        for i in range(n):
        	indice=np.int32((e[i]-emin)/(emax-emin)*(nbines-1))
        	if( (jz[i] >= yy[indice]) & (rho[i] >= threshold) ):
        		yy[indice]=jz[i] 
        
        k=np.where(yy > -1e-10)
        n=np.size(k)
        #n=n[0]
        
        ecirb=xx[k]
        jcirb=yy[k]
        
        #print min(e1)
        
        
        print '-------------------------------------------------'
        
        ejcirb=np.ndarray([np.size(ecirb),2])
        ejcirb[:,0]=ecirb
        ejcirb[:,1]=jcirb
        
        jcir=np.ndarray([np.size(e)])
        
        k=1.
        for j in range(np.size(e)):
                #idx = np.argmin(np.abs(e[j] - ecirb))
                kk=np.where((e[j]<ecirb))
                #print kk
                if np.size(kk)>0:
                        kk0=min(kk[0])
                                
                else:
                        kk0=np.size(ecirb)-1
        		jcir[j]=max(jcirb)*1.01
        		continue
        
        
        	if kk0==0:
                        kk0=kk0+1
        
                kk1=kk0
                kk0=kk1-1
                a=(jcirb[kk1]-jcirb[kk0])/(ecirb[kk1]-ecirb[kk0])
                b=(((jcirb[kk0]*ecirb[kk1])-(jcirb[kk1]*ecirb[kk0]))/(ecirb[kk1]-ecirb[kk0]))
                jcir[j]=e[j]*a+b
        

	return jcir
