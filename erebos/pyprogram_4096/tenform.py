import numpy as np

def tenf(x,y,z):
        size=np.size(x)
        a00=np.sum(x*x)
        a01=np.sum(x*y)
        a02=np.sum(x*z)
        a10=np.sum(y*x)
        a11=np.sum(y*y)
        a12=np.sum(y*z)
        a20=np.sum(z*x)
        a21=np.sum(z*y)
        a22=np.sum(z*z)

        T=np.matrix([[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]])
        T=T/size                 
        return T

def tenfr(x,y,z,m):

        r2=(x**2+y**2+z**2)
        size=sum(m)
        a00=np.sum(m*x*x/r2)
        a01=np.sum(m*x*y/r2)
        a02=np.sum(m*x*z/r2)
        a10=np.sum(m*y*x/r2)
        a11=np.sum(m*y*y/r2)
        a12=np.sum(m*y*z/r2)
        a20=np.sum(m*z*x/r2)
        a21=np.sum(m*z*y/r2)
        a22=np.sum(m*z*z/r2)

        T=np.matrix([[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]])
        T=T/size                 
        return T

def tenf2d(x,y):
        size=np.size(x)
        a00=np.sum(x*x)
        a01=np.sum(x*y)
        a11=np.sum(y*y)
        a10=np.sum(y*x)
        T=np.matrix([[a00, a01], [a10, a11]])
        T=T/size                 
        return T
    
def tenfr2d(x,y):
        size=np.size(x)
        r2=x**2+y**2
        a00=np.sum(x*x/r2)
        a01=np.sum(x*y/r2)
        a11=np.sum(y*y/r2)
        a10=np.sum(y*x/r2)
        T=np.matrix([[a00, a01], [a10, a11]])
        T=T/size                 
        return T
