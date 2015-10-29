import numpy as np
#Some constants
G = 4.52e-48 #Newton's gravitional constant in Mpc^3/s^2/Solar Mass
Mpcperkm = 3.241e-20 #Mpc/km; used to convert H0 to s^-1

def get_nfw_profile(r_s,z,M,r,H0,OM,delta):
    rhocrit = 3*(H0*Mpcperkm)**2/(8*np.pi*G)/(H0/100.)**2
    #rhocrit has units of SM h^2/Mpc^3
    rhom = rhocrit*OM
    rdelta = (M/(4./3.*np.pi*rhom*delta))**(1./3.)#The virial radius in Mpc/h
    c = rdelta/r_s#The concentration
    fc = np.log(1+c)-c/(1+c)
    rrs = r/r_s #the radii divided by the scale radius
    profile = np.zeros((len(z),len(r)))
    #Loop over z is in the :
    profile[:] = M/(4.*np.pi*r_s**3*fc)/(rrs*(1+rrs)**2)
    #NOTE!!!: HW08 have a factor of OM in the denominator but Eduardo doesn't
    return [profile,c] #Units of SM h^2/Mpc^3
