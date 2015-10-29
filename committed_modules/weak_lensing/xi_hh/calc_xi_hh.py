import numpy as np

#Multiply the bias onto xi_mm
def calc_xi_hh(bias,xi_mm):
    xi_hh = np.zeros(np.shape(xi_mm))
    for i in xrange(0,len(bias)):
        xi_hh[i] = bias[i]*bias[i]*xi_mm[i]
    return xi_hh
    
