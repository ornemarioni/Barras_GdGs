def set_tick(ax, size=16):
    ax.minorticks_on()
    ax.tick_params(labelsize=size)
    ax.tick_params('both', length=3, width=1.2,which='minor', direction='in', right='on',top='on')
    ax.tick_params('both', length=6, width=1.2,which='major', direction='in', right='on',top='on') 
    ax.tick_params('both', length=6, width=1.2,which='major', direction='in', right='on',top='on')
    return ax