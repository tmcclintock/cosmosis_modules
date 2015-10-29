import numpy as np

#Multiply the bias onto xi_mm
def calc_xi_2halo(bias,xi_mm):
    xi_2halo = np.zeros(np.shape(xi_mm))
    for i in xrange(0,len(bias)):
        xi_2halo[i] = bias[i]*xi_mm[i]
    return xi_2halo
    
