import sphviewer as sph


##posiciones de partículas que se quiere graficar

pos=np.ndarray([3,np.size(x4n)])
pos[0,:]=x4n
pos[1,:]=y4n
pos[2,:]=z4n

#generador del gráfico
particles=sph.Particles(pos,m4,nb=5)

#escena del gráfico---------
escena=sph.Scene(particles)

#como queres la escena-----------
escena.set_autocamera(mode='density')
rl=20
escena.update_camera(r='infinity',x=0,y=0,z=0,extent=[-rl,rl,-rl,rl])

#renderizado de la escena (aca hace la grid y cuenta particulas)
rend=sph.Render(escena)

##extencion del grafico
extent=escena.get_extent()

#escala logaritmica-----------------
rend.set_logscale()

#rango que tiene la escala  de colores-----
vmin=4
vmax=9

# escala de colores que te guste (http://matplotlib.org/examples/color/colormaps_reference.html)
cmap='jet'

#grafico-----------
fig=plt.figure(5,figsize=(10,10))
ax=fig.add_subplot(111)
ax.plot(rend.get_image(),extent=extent,vmin=vmin,vmax=vmax,origin='lower',cmap='jet')

#etc (como queres el grafico depende de vos)
