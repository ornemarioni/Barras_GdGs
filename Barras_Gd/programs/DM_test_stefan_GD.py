import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import rotation_mio as rot
import numpy as np
import barstrength2 as strng
import time_conversion as tiempo
import sphviewer as sph
# from mpl_toolkits.axes_grid1 import make_axes_locatable


h=0.73005
G = 4.299e-6
a0=1
H0 = h*100
omega_lambda=0.716
omega_matter=0.1277/(h**2.)
omega0 = omega_lambda + omega_matter
vector2 = ('M31','MW','M33')
vector = ('00','01','03')
vector3 = ('A','B')
# carpeta = ('9in1_M31/','9in1_MW_new/')
# carpeta = ('test_stefan/')


#path = 'home/ornela/SimCLUES/'
path = '/home/omarioni/'

atime = np.loadtxt(path + 'redshift_outputs.txt')
aexp = atime[:,2]
# aexp = aexp[::-1]

path = '/mnt/is2/alejandro/ornella/'

snapshot=range(496,100,-1)
for isnap in snapshot:
    snap = h5py.File(path + 'outputs_1/snap_'+str('%03d'%isnap)+'.h5py', 'r')
    print isnap
#     for i in range(1,2):

    i=1
    
    cm   = snap['subhalo_0'+ str('%s' %vector[i])+ '/Center'][()]
    r200 = snap['subhalo_0'+ str('%s' %vector[i])+ '/R200'][()]

    pstr = snap['subhalo_0'+ str('%s' %vector[i]) + '/Str/Coordinates'][()]
    mstr = snap['subhalo_0'+ str('%s' %vector[i]) + '/Str/Masses'][()]
    vel  = snap['subhalo_0'+ str('%s' %vector[i])+ '/Str/Velocities'][()]

    pdrk = snap['subhalo_0'+ str('%s' %vector[i]) + '/Drk/Coordinates'][()]
    mdrk = snap['subhalo_0'+ str('%s' %vector[i]) + '/Drk/Masses'][()]

    z = a0/aexp[isnap] - 1.
    Ht = H0*np.sqrt(omega_lambda+(1-omega0)*(1+z)**2+omega_matter*(1+z)**3)
    time = tiempo.conv(z, h, omega_lambda, omega_matter)

#         print isnap#, time

    #---aca paso las coordenadas respecto al centro de la galaxia------
    xstr = (pstr[:,0]-cm[0])*aexp[isnap]/h
    ystr = (pstr[:,1]-cm[1])*aexp[isnap]/h
    zstr = (pstr[:,2]-cm[2])*aexp[isnap]/h
    r = np.sqrt(xstr**2+ystr**2+zstr**2)

    xdrk = (pdrk[:,0]-cm[0])*aexp[isnap]/h
    ydrk = (pdrk[:,1]-cm[1])*aexp[isnap]/h
    zdrk = (pdrk[:,2]-cm[2])*aexp[isnap]/h
    rdrk = np.sqrt(xdrk**2+ydrk**2+zdrk**2)

    v_x = vel[:,0] *np.sqrt(aexp[isnap]) + Ht *xstr/1000.
    v_y = vel[:,1] *np.sqrt(aexp[isnap]) + Ht *ystr/1000.
    v_z = vel[:,2] *np.sqrt(aexp[isnap]) + Ht *zstr/1000.

    #----------------------masas----------------------------
    mstr = mstr/h
    mdrk = mdrk/h

    r200 = r200*aexp[isnap]/h
    rgal=0.15*r200

    limit = np.where(r<rgal)
    r_sort = np.sort(r[limit])
    r_indice = np.argsort(r[limit])

    Mc_str = np.cumsum((mstr[limit])[r_indice])
    M_gal = Mc_str[-1]

    #------------------ calculamos M90------------------------------------------
    razon = Mc_str/M_gal
    noventa, = np.where(razon < 0.9)
    cincuenta, = np.where(razon < 0.5)

    r90 = r_sort[noventa][-1]
    r50 = r_sort[cincuenta][-1]             

    #--------------------------------------------         
    veloc,=np.where(r<r50)

    #----------componentes de la velocidad del centro de masa------------
    vxcm = sum(mstr[veloc]*v_x[veloc])/sum(mstr[veloc])
    vycm = sum(mstr[veloc]*v_y[veloc])/sum(mstr[veloc])
    vzcm = sum(mstr[veloc]*v_z[veloc])/sum(mstr[veloc])

    #----- velocidades de las estrellas respecto del centro de masa de la galaxia---------
    vx = v_x - vxcm
    vy = v_y - vycm
    vz = v_z - vzcm

    e1x,e2x,e3x,e1y,e2y,e3y,e1z,e2z,e3z = rot.rot1(mstr,xstr,ystr,zstr,vx,vy,vz,3*aexp[isnap])

# posiciones de particulas que se quiere graficar
# como lo de arriba me da los versores hago las posiciones con esto 

    xn = e1x*xstr + e1y*ystr + e1z*zstr
    yn = e2x*xstr + e2y*ystr + e2z*zstr
    zn = e3x*xstr + e3y*ystr + e3z*zstr
    vxn = e1x*vx + e1y*vy + e1z*vz
    vyn = e2x*vx + e2y*vy + e2z*vz
    vzn = e3x*vx + e3y*vy + e3z*vz

    xn_drk = e1x*xdrk + e1y*ydrk + e1z*zdrk
    yn_drk = e2x*xdrk + e2y*ydrk + e2z*zdrk
    zn_drk = e3x*xdrk + e3y*ydrk + e3z*zdrk


    pos=np.ndarray([np.size(xn),4])
    pos[:,0]=xn
    pos[:,1]=yn
    pos[:,2]=zn
    pos[:,3]=mstr

    pos2=np.ndarray([np.size(xn_drk),4])
    pos2[:,0]=xn_drk
    pos2[:,1]=yn_drk
    pos2[:,2]=zn_drk
    pos2[:,3]=mdrk


    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(16, 16))

    fig.subplots_adjust(bottom=0.07, left =0.09, right = 0.97, top = 0.95, wspace=0.0, hspace= 0.0)

    #----------------------------------------------------------------------

    #---------------------generador del grafico3-----------------
    rl= 100   
    corte,=np.where((xn_drk <rl) & (yn_drk <rl) & (zn_drk <rl) & (xn_drk >-rl) & (yn_drk >-rl) & (zn_drk >-rl))


    #-----rango que tiene la escala  de colores-----
    vmin=4.5
    vmax=6.7

    # ----escala de colores que te guste (http://matplotlib.org/examples/color/colormaps_reference.html)---
    cmap='viridis'

    #         nb1 = 300 
    nb1 = 25
    npixel = 1000

    particles=sph.Particles(pos2[corte,:3],mdrk[corte]*1e10,nb=nb1)
    escena=sph.Scene(particles)
    escena.update_camera(r='infinity',x=0,y=0,z=0,extent=[-rl,rl,-rl,rl], xsize=npixel,ysize=npixel)
    rend=sph.Render(escena)
    extent=escena.get_extent()
    rend.set_logscale()

    ax[1,0].imshow(rend.get_image(),extent=extent,origin='lower',cmap=cmap, vmin=vmin, vmax= vmax)
    ax[1,0].set_xlim(-100,100)
    ax[1,0].set_ylim(-100,100)
    ax[1,0].set_xticks([-100,-50,0,50,100])
    ax[1,0].set_yticks([-100,-50,0,50,100])
    ax[1,0].set_xlabel('x',fontsize=28)
    ax[1,0].set_ylabel('y',fontsize=28)
    ax[1,0].minorticks_on()
    ax[1,0].tick_params( labelsize=24)
    ax[1,0].tick_params('both', length=10, width=2,which='minor', direction='in', right=True,top=True)
    ax[1,0].tick_params('both', length=15, width=2,which='major', direction='in', right=True,top=True)
    ax[1,0].add_patch(plt.Circle((0,0),radius=rgal,ec='k',fc=None,ls='--',lw=2.5, fill=False))
    #     ax[0,1].set_xticklabels([])
    #     ax[0,1].set_yticklabels([])
    #     ax[0,1].annotate("",xy=(-40, -40), xycoords='data',xytext=(-10, -40),textcoords='data',
    #                  ha='center', va='center', 
    #                 arrowprops=dict(arrowstyle="|-|", connectionstyle='arc3', color ='white', lw=2.5))

    #     ax[2,0].text(-25, -40, '30kpc', fontsize=25, color='white', ha='center', va='bottom')


    #--------------------------------------
    particles=sph.Particles(pos2[corte,:3],mdrk[corte]*1e10,nb=nb1)
    escena=sph.Scene(particles)
    escena.update_camera(r='infinity',x=0,y=0,z=0,extent=[-rl,rl,-rl,rl], t=90, xsize=npixel,ysize=npixel)
    rend=sph.Render(escena)
    extent=escena.get_extent()
    rend.set_logscale()


    ax[0,0].imshow(rend.get_image(),extent=extent,origin='lower',cmap=cmap, vmin=vmin, vmax= vmax)
    ax[0,0].set_xlim(-100,100)
    ax[0,0].set_ylim(-100,100)
    ax[0,0].set_xticks([-100,-50,0,50,100])
    ax[0,0].set_yticks([-50,0,50,100])
    ax[0,0].set_xticklabels([])
    ax[0,0].set_yticklabels([-50,0,50,100])
    ax[0,0].set_ylabel('z',fontsize=28)
    ax[0,0].minorticks_on()
    ax[0,0].tick_params( labelsize=24)
    ax[0,0].tick_params('both', length=10, width=2,which='minor', direction='in', right=True,top=True)
    ax[0,0].tick_params('both', length=15, width=2,which='major', direction='in', right=True,top=True)
    ax[0,0].text(-85, 82,'z='+str('%.3f'%z), fontsize=25, color='yellow', ha='left', va='center') 
    ax[0,0].add_patch(plt.Circle((0,0),radius=rgal,ec='k',fc=None,ls='--',lw=2.5, fill=False))
    ax[0,0].set_title('B-GADGET',fontsize=25,loc='center')

#     #--------------------------------------
    particles=sph.Particles(pos2[corte,:3],mdrk[corte]*1e10,nb=nb1)
    escena=sph.Scene(particles)
    escena.update_camera(r='infinity',x=0,y=0,z=0,extent=[-rl,rl,-rl,rl], t=90,p=90, xsize=npixel,ysize=npixel)
    rend=sph.Render(escena)
    extent=escena.get_extent()
    rend.set_logscale()


    ax[0,1].imshow(rend.get_image(),extent=extent,origin='lower',cmap=cmap, vmin=vmin, vmax= vmax)
    ax[0,1].set_xlim(-100,100)
    ax[0,1].set_ylim(-100,100)
    ax[0,1].set_xticks([-50,0,50,100])
    ax[0,1].set_yticks([-100,-50,0,50,100])
    ax[0,1].set_yticklabels([])
    ax[0,1].set_xticklabels([-50,0,50,100])
    ax[0,1].minorticks_on()
    ax[0,1].tick_params( labelsize=24)
    ax[0,1].tick_params('both', length=10, width=2,which='minor', direction='in', right=True,top=True)
    ax[0,1].tick_params('both', length=15, width=2,which='major', direction='in', right=True,top=True)
    ax[0,1].set_xlabel('y',fontsize=28)
    ax[0,1].text(85, 82, str('%.3f'%time)+'Gyr', fontsize=25, color='yellow', ha='right', va='center') 
    ax[0,1].add_patch(plt.Circle((0,0),radius=rgal,ec='k',fc=None,ls='--',lw=2.5, fill=False))


    ax[1,1].axis('off')
    
    plt.show(False)
    
    path2 = '/home/omarioni/Barras_GdGs/Barras_Gd/_imagenes/test_STEFAN/'
    fig.savefig(path2 + str('%s' %vector2[i])+'_'+str('%s' %isnap)+'.png',
                dpi = 100, xxbox_inches='tight')
    plt.close()
