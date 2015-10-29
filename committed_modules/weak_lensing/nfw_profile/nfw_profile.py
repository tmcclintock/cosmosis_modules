from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
import numpy as np
import calc_nfw

cosmo = section_names.cosmological_parameters
distances = section_names.distances
def setup(options):
    #Calculate the radii at which we will calculate the correlation function
    #Note: we need to add a slight offset to the upper bound
    #on r so that we are gauranteed to get numr radii
    rmin = options[option_section,"r_min"]
    rmax = options[option_section,"r_max"]
    numr = options[option_section,"num_r"]
    dlr = (np.log(rmax)-np.log(rmin))/(numr-1)
    r = np.exp((np.arange(np.log(rmin),np.log(rmax)+0.00001,dlr)))

    #Get the overdensity parameter
    delta = options[option_section,"delta"]

    return (delta,r)

def execute(block,config):
    #Take in the configured parameters
    delta = config[0]
    r = config[1]

    #Acquire the scale radius and the cluster mass
    r_s = block[cosmo,"scale_radius"]
    M = block[cosmo,"log10_cluster_mass"]
    M = np.power(10,M)

    #Acquire the redshift samples
    z = block[distances,"z"]

    #Take in the cosmological parameters
    H0 = block[cosmo,"hubble"]
    OM = block[cosmo,"omega_m"]

    #Calculate the nfw profile
    nfw_profile,c = calc_nfw.get_nfw_profile(r_s,z,M,r,H0,OM,delta)

    #Save to the block
    section = "nfw_profile"
    block[section,"concentration"]=c
    block[section,"rho_nfw"]=nfw_profile

    return 0
        
def cleanup(self):
    #No cleanup necessary
    return 0


