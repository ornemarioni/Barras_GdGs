import h5py
snap = h5py.File('/mnt/is2/alejandro/ornella/outputs/snap_497.h5py', 'r')
def printname(name):
    print name
snap.visit(printname)

