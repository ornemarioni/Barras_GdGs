import numpy as np
import bines
import math



def a2max(m,x,y,nbin):


    rcil=np.sqrt(x**2+y**2)
    rbin=bines.rbin1(rcil,np.ones(np.size(m)),nbin)
    A2v=np.ndarray(nbin)
    a2=np.ndarray(nbin)
    b2=np.ndarray(nbin)
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
        a2[j]=sum(m[kk]*np.cos(2*titaj))
        b2[j]=sum(m[kk]*np.sin(2*titaj))
        A2v[j]=np.sqrt(a2[j]**2+b2[j]**2)/a0

    A2max=max(A2v)
    kA2=np.argmax(A2v) 
    rbinmax=rbin[kA2]
    a2_max=a2[kA2]
    b2_max=b2[kA2]

    return a2_max, b2_max, rbinmax

def a2(m,x,y,nbin):

    rcil=np.sqrt(x**2+y**2)
    rbin=bines.rbin1(rcil,np.ones(np.size(m)),nbin)
    A2v=np.ndarray(nbin)
    phiv=np.ndarray(nbin)
    a2=np.ndarray(nbin)
    b2=np.ndarray(nbin)
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
        a2[j]=sum(m[kk]*np.cos(2*titaj))
        b2[j]=sum(m[kk]*np.sin(2*titaj))
        A2v[j]=np.sqrt(a2[j]**2+b2[j]**2)/a0


    return A2v,rbin, a2, b2
