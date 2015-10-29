import numpy as np

#Some constants
G = 4.52e-48 #Newton's gravitional constant in Mpc^3/s^2/Solar Mass
Mpcperkm = 3.241e-20 #Mpc/km; used to convert H0 to s^-1

#Calculate the 1 halo correlation function
def calc_xi_1halo(rho,H0,OM):
    #rho has units of SM h^2/Mpc^3
    rhocrit = 3*(H0*Mpcperkm)**2/(8*np.pi*G)/(H0/100.)**2
    #rhocrit has units of SM h^2/Mpc^3
    rhom = OM*rhocrit#Units of SM h^2/Mpc^3
    xi_1halo = rho/rhom - 1
    return xi_1halo
