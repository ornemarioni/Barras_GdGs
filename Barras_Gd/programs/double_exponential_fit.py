import numpy as np
from scipy.optimize import curve_fit

'''This program calculate a double exponential fit for a dataset'''

def func(x, a, b, c, d):
    f = np.log(a * np.exp(b*x) + c * np.exp(d*x))
    return f

def double_exp_fit(x, y, cut_1, cut_2):
    
    #This fit 2 linear polynomials to generate the initial params of the
    #double exponential fit
    
    cut1, = np.where((x > cut_1[0]) & (x < cut_1[1]))
    linear1 = np.polyfit(x[cut1], np.log(y[cut1]), 1)
    
    cut2, = np.where((x > cut_2[0]) & (x < cut_2[1]))
    linear2 = np.polyfit(x[cut2], np.log(y[cut2]), 1)
    
    
    param = np.asarray([np.exp(linear1[1]), -1/linear1[0], 
                        np.exp(linear2[1]), -1/linear2[0]])
    
    cut, = np.where((x > cut_1[0]) & (x < cut_2[1]))
    
    popt, pcov = curve_fit(func, x[cut], np.log(y[cut]), p0=param)
    
    return popt, pcov