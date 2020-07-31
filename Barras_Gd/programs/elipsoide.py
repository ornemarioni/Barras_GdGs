## Calculamos las particulas dentro de un elipsoide de semieje mayor = radio
## el angulo es el phi, es para alinear el semieje mayor con el eje x
import numpy as np
import tenform as ten
import bines2 as bine
import scipy.interpolate as sint

def elipsoide(x,y,z,vx,vy,vz,radio,phi=0, nbin=20):
    
    if radio == np.nan:
        omegabar = np.nan
    else:
        xx = x* np.cos(phi) + y*np.sin(phi)
        yy = x*-np.sin(phi) + y*np.cos(phi)
        zz = z

        r = np.sqrt(xx**2 + yy**2 + zz**2)

        limit, = np.where(r < radio)

        tensor   = ten.tenf(xx[limit], yy[limit], zz[limit]) #calculo el tensor de forma
        matriz   = np.linalg.eig(tensor) #saco los autovalores
        autoval  = matriz[0]
        autovec  = matriz[1]
        
#         print autovec, '\n'

        asort  = np.sort(autoval) #los ordeno de menor a mayor

        A = np.sqrt(asort[2]) #semiejes
        B = np.sqrt(asort[1])
        C = np.sqrt(asort[0])

        aa = radio #normalizo al radio que yo quiero
        bb = (B/A)*radio
        cc = (C/A)*radio

        Relip = np.sqrt((xx/aa)**2 + (yy/bb)**2 + (zz/cc)**2) # formula del elipsoide
        mask, = np.where(Relip <= 1.)

        xn  = x[mask]
        yn  = y[mask]
        zn  = z[mask]
        vxn = vx[mask]
        vyn = vy[mask]
        vzn = vz[mask]

        rn = np.sqrt(xn**2 + yn**2 + zn**2)
        Rn = np.sqrt(xn**2 + yn**2)

        Vtg = (-yn*vxn + xn*vyn)/Rn

        omega = Vtg/Rn

        med, nodos = bine.rbin1(rn,nbin)
        omega_mean = np.zeros(nbin)

        for i in range(0,nbin):
            inbin, = np.where((rn > nodos[i]) & (rn < nodos[i+1]))

            omega_mean[i] = np.mean(omega[inbin])

        finterp = sint.interp1d(med,omega_mean,fill_value="extrapolate")

        omegabar = finterp(radio)
    
    return omegabar, xn, yn, zn, autovec