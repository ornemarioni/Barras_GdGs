import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import rotation_mio as rot
import numpy as np
import barstrength2 as strng
import time_conversion as tiempo
import sphviewer as sph
from mpl_toolkits.axes_grid1 import make_axes_locatable

a0 = 1.
vector = (1,2,4)
vector2 = ('M31', 'MW')
vector3 = ('A','B')
# carpeta = ('9in1_M31/', '9in1_MW_new/')

#path = 'home/ornela/SimCLUES/'
path = '/mnt/sersic2/omarioni/'

snapshot = np.loadtxt(path + 'Gasoline/snapshots.txt', dtype='string')


for j in range(59,0,-1):
    snap = h5py.File(path + 'Gasoline/outputs2/snap_'+str('%s'%snapshot[j])+'.h5py', 'r')

    print snapshot[j]

#     for i in range(1,2):
    i=1
    
    cm   = snap['subhalo_00'+ str('%s' %vector[i])+ '/Center'].value
    r200 = snap['subhalo_00'+ str('%s' %vector[i])+ '/R200'].value
    time = snap['subhalo_00'+ str('%s' %vector[i])+ '/Time'].value
    h    = snap['subhalo_00'+ str('%s' %vector[i])+ '/h'].value
    aexp = snap['subhalo_00'+ str('%s' %vector[i])+ '/aexp'].value

    pstr = snap['subhalo_00'+ str('%s' %vector[i]) + '/Str/Coordinates'].value
    mstr = snap['subhalo_00'+ str('%s' %vector[i]) + '/Str/Masses'].value
    vel  = snap['subhalo_00'+ str('%s' %vector[i])+ '/Str/Velocities'].value

    pdrk = snap['subhalo_00' + str('%s' %vector[i]) + '/Drk/Coordinates'].value
    mdrk = snap['subhalo_00' + str('%s' %vector[i]) + '/Drk/Masses'].value

    z = a0/aexp - 1.
    #---aca paso las coordenadas respecto al centro de la galaxia------
    xstr = pstr[:,0]-cm[0]
    ystr = pstr[:,1]-cm[1]
    zstr = pstr[:,2]-cm[2]
    r = np.sqrt(xstr**2+ystr**2+zstr**2)

    xdrk = pdrk[:,0]-cm[0]
    ydrk = pdrk[:,1]-cm[1]
    zdrk = pdrk[:,2]-cm[2]
    rdrk = np.sqrt(xdrk**2+ydrk**2+zdrk**2)

    v_x = vel[:,0] 
    v_y = vel[:,1] 
    v_z = vel[:,2] 

    #----------------------masas----------------------------

    r200 = r200*aexp
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

    if j==59:
        e1x,e2x,e3x,e1y,e2y,e3y,e1z,e2z,e3z = rot.rot1(mstr,xstr,ystr,zstr,vx,vy,vz,3*aexp)

    ##posiciones de particulas que se quiere graficar
    ##como lo de arriba me da los versores hago las posiciones con esto 

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
    rl= 50   
    corte,=np.where((xn <rl) & (yn <rl) & (zn <rl) & (xn >-rl) & (yn >-rl) & (zn >-rl))


    #-----rango que tiene la escala  de colores-----
    vmin=0.9
    vmax=6.7

    # ----escala de colores que te guste (http://matplotlib.org/examples/color/colormaps_reference.html)---
    cmap='magma'

    #         nb1 = 300 
    nb1 = 25
    npixel = 1000

    particles=sph.Particles(pos[corte,:3],mstr[corte],nb=nb1)
    escena=sph.Scene(particles)
    escena.update_camera(r='infinity',x=0,y=0,z=0,extent=[-rl,rl,-rl,rl], xsize=npixel,ysize=npixel)
    rend=sph.Render(escena)
    extent=escena.get_extent()
    rend.set_logscale()

    ax[1,0].imshow(rend.get_image(),extent=extent,origin='lower',cmap=cmap, vmin=vmin, vmax= vmax)
    ax[1,0].set_xlim(-50,50)
    ax[1,0].set_ylim(-50,50)
    ax[1,0].set_xticks([-50,-25,0,25,50])
    ax[1,0].set_yticks([-50,-25,0,25,50])
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
    particles=sph.Particles(pos[corte,:3],mstr[corte],nb=nb1)
    escena=sph.Scene(particles)
    escena.update_camera(r='infinity',x=0,y=0,z=0,extent=[-rl,rl,-rl,rl], t=90, xsize=npixel,ysize=npixel)
    rend=sph.Render(escena)
    extent=escena.get_extent()
    rend.set_logscale()


    ax[0,0].imshow(rend.get_image(),extent=extent,origin='lower',cmap=cmap, vmin=vmin, vmax= vmax)
    ax[0,0].set_xlim(-50,50)
    ax[0,0].set_ylim(-50,50)
    ax[0,0].set_xticks([-50,-25,0,25,50])
    ax[0,0].set_yticks([-25,0,25,50])
    ax[0,0].set_xticklabels([])
    ax[0,0].set_yticklabels([-25,0,25,50])
    ax[0,0].set_ylabel('z',fontsize=28)
    ax[0,0].minorticks_on()
    ax[0,0].tick_params( labelsize=24)
    ax[0,0].tick_params('both', length=10, width=2,which='minor', direction='in', right=True,top=True)
    ax[0,0].tick_params('both', length=15, width=2,which='major', direction='in', right=True,top=True)
    ax[0,0].text(-45, 45,'z='+str('%.3f'%z), fontsize=25, color='yellow', ha='left', va='center') 
    ax[0,0].add_patch(plt.Circle((0,0),radius=rgal,ec='k',fc=None,ls='--',lw=2.5, fill=False))
    ax[0,0].set_title('B-GASOLINE',fontsize=25,loc='center')

#     #--------------------------------------
    particles=sph.Particles(pos[corte,:3],mstr[corte],nb=nb1)
    escena=sph.Scene(particles)
    escena.update_camera(r='infinity',x=0,y=0,z=0,extent=[-rl,rl,-rl,rl], t=90,p=90, xsize=npixel,ysize=npixel)
    rend=sph.Render(escena)
    extent=escena.get_extent()
    rend.set_logscale()


    ax[0,1].imshow(rend.get_image(),extent=extent,origin='lower',cmap=cmap, vmin=vmin, vmax= vmax)
    ax[0,1].set_xlim(-50,50)
    ax[0,1].set_ylim(-50,50)
    ax[0,1].set_xticks([-25,0,25,50])
    ax[0,1].set_yticks([-50,-25,0,25,50])
    ax[0,1].set_yticklabels([])
    ax[0,1].set_xticklabels([-25,0,25,50])
    ax[0,1].minorticks_on()
    ax[0,1].tick_params( labelsize=24)
    ax[0,1].tick_params('both', length=10, width=2,which='minor', direction='in', right=True,top=True)
    ax[0,1].tick_params('both', length=15, width=2,which='major', direction='in', right=True,top=True)
    ax[0,1].set_xlabel('y',fontsize=28)
    ax[0,1].text(45, 45, str('%.3f'%time)+'Gyr', fontsize=25, color='yellow', ha='right', va='center') 
    ax[0,1].add_patch(plt.Circle((0,0),radius=rgal,ec='k',fc=None,ls='--',lw=2.5, fill=False))


    ax[1,1].axis('off')
    
    plt.show(False)
    
    path2 = '/home/omarioni/Barras_GdGs/Barras_Gs/_imagenes/test_STEFAN_50kpc/STAR/'
    fig.savefig(path2 + str('%s' %vector2[i])+'_'+str('%s' %snapshot[j])+'.png',
                dpi = 100, xxbox_inches='tight')
    plt.close()