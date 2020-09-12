import numpy as np

'''This program calculate the instability criterium of '''
def f_disk(velocity, mass, radius):
    f_disk = velocity / np.sqrt(G*mass/radius)