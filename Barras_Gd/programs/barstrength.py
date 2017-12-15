import numpy as np
import bines
import math



def a2max(m,x,y,nbin):


    rcil=np.sqrt(x**2+y**2)
    rbin=bines.rbin1(rcil,np.ones(np.size(m)),nbin)
    A2v=np.ndarray(nbin)
    phiv=np.ndarray(nbin)
    out=np.ndarray([2,nbin])
    for j in range(1,nbin):
        kk,=np.where((rcil>rbin[j-1]) & (rcil<=rbin[j]))
        xnk=x[kk]
        ynk=y[kk]
        titaj=np.ndarray(np.size(kk))
        for ii in range(np.size(kk)):
            titaj[ii]=math.atan2(ynk[ii],xnk[ii])

        rj=np.sqrt(x[kk]**2+y[kk]**2) 
        a0=sum(m[kk])
        a2=sum(m[kk]*np.cos(2*titaj))
        b2=sum(m[kk]*np.sin(2*titaj))
        A2v[j]=np.sqrt(a2**2+b2**2)/a0
        phiv[j]=math.atan2(b2,a2)/2.

    A2max=max(A2v)
    kA2=np.argmax(A2v) 
    rbinmax=rbin[kA2]
    phimax=phiv[kA2]

    return A2max,rbinmax, phimax

def a2(m,x,y,nbin):

    rcil=np.sqrt(x**2+y**2)
    rbin=bines.rbin1(rcil,np.ones(np.size(m)),nbin)
    A2v=np.ndarray(nbin)
    phiv=np.ndarray(nbin)
    out=np.ndarray([2,nbin])
    for j in range(1,nbin):
        kk,=np.where((rcil>rbin[j-1]) & (rcil<=rbin[j]))
        xnk=x[kk]
        ynk=y[kk]
        titaj=np.ndarray(np.size(kk))
        for ii in range(np.size(kk)):
            titaj[ii]=math.atan2(ynk[ii],xnk[ii])

        rj=np.sqrt(x[kk]**2+y[kk]**2) 
        a0=sum(m[kk])
        a2=sum(m[kk]*np.cos(2*titaj))
        b2=sum(m[kk]*np.sin(2*titaj))
        A2v[j]=np.sqrt(a2**2+b2**2)/a0
        phiv[j]=math.atan2(b2,a2)/2.

    A2max=max(A2v)
    kA2=np.argmax(A2v) 
    rbinmax=rbin[kA2]

    return A2v,phiv,rbin
