from __future__ import unicode_literals

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True

fig=plt.figure(1, figsize=(8,8))
fig.subplots_adjust(bottom=0.12, left =0.16, right = 0.95, top = 0.95)
ax=fig.add_subplot(111)

ax.plot(lbar2_M31gd[0],rcor_phi[0],'*', markersize=25, mec='b',mfc='b', mew=3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.25,70)
ax.set_ylim(0.25,70)
ax.minorticks_on()
ax.tick_params( labelsize=22)
ax.tick_params('both', length=3, width=1.2,which='minor', direction='in', right='on',top='on')#,colors='w')
ax.tick_params('both', length=6, width=1.2,which='major', direction='in', right='on',top='on')#,colors='w')  
ax.set_xlabel(r'$l_{bar}\:[kpc]$', fontsize=30)#, color='lightgray')
ax.set_ylabel(r'$R_{corot}\:[kpc]$',fontsize=30)#, color='lightgray')
ax.text(x = 0.4, y = 30, s = u'Slow Bars', fontsize = 30, va = 'bottom', ha = 'left', color='k', style='italic')
ax.text(x = 15, y = 4, s = u'Fast Bars', fontsize = 30, va = 'bottom', ha = 'left',  color='k', style='italic')
#ax.legend([algorry, corsini, font],['Algorry et al. 2017', 'Corsini 2011', 'Font et al. 2017'],fontsize=25, frameon=False, loc=4)
ax.arrow( lbar2_M31gd[109], rcor_phi[109], lbar2_M31gd[0]-lbar2_M31gd[109] , rcor_phi[0]-rcor_phi[109],
         fc='k', ec='k',lw=5, head_width=0.2, head_length=0.3, zorder=10 )

fig.savefig('/home/omarioni/Barras_GdGs/Barras_Gd/_imagenes/Rcor_lbar.pdf', dpi = 100, xxbox_inches='tight')
plt.show()