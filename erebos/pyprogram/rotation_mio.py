import numpy as np

def rot1(m,x,y,z,vx,vy,vz,corte):

        r   = np.sqrt(x**2+y**2+z**2)
        kk, = np.where(r<corte)        
    
        jx = sum(m[kk]*(y[kk]*vz[kk]-z[kk]*vy[kk]))
        jy = sum(m[kk]*(z[kk]*vx[kk]-x[kk]*vz[kk]))
        jz = sum(m[kk]*(x[kk]*vy[kk]-y[kk]*vx[kk]))

        j = np.sqrt(jx**2+jy**2+jz**2)

        rjx = jx/j
        rjy = jy/j
        rjz = jz/j

        jp = np.sqrt(jx**2+jy**2)

        e1x = jy/jp
        e1y = -jx/jp
        e1z = 0.

        #e1=np.sqrt(e1x**2+e1y**2+e1z**2)
        #e1x=e1x/e1
        #e1y=e1y/e1
        #e1z=e1z/e1

        e2x = jx*jz/(j*jp)
        e2y = jy*jz/(j*jp)
        e2z = -jp/j
        #e2=np.sqrt(e2x**2+e2y**2+e2z**2)
        #e2x=e2x/e2
        #e2y=e2y/e2
        #e2z=e2z/e2

        e3x = rjx
        e3y = rjy
        e3z = rjz
        #e3=np.sqrt(e3x**2+e3y**2+e3z**2)
        #e3x=e3x/e3
        #e3y=e3y/e3
        #e3z=e3z/e3

        xn=e1x*x+e1y*y+e1z*z
        yn=e2x*x+e2y*y+e2z*z
        zn=e3x*x+e3y*y+e3z*z
        vxn=e1x*vx+e1y*vy+e1z*vz
        vyn=e2x*vx+e2y*vy+e2z*vz
        vzn=e3x*vx+e3y*vy+e3z*vz

        return e1x,e2x,e3x,e1y,e2y,e3y,e1z,e2z,e3z 
