from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
import numpy as np
import calc_xi_1halo

cosmo = section_names.cosmological_parameters
distances = section_names.distances
def setup(options):
    #Take in the profile used
    prof = options[option_section,"profile_name"]
    return (prof)

def execute(block,config):
    #Take in the configured parameters
    prof = config
    profile_name = prof+"_profile"
    #this is the actual name of the profile blocks
    
    #Take in the cosmological parameters
    H0 = block[cosmo,"hubble"]
    OM = block[cosmo,"omega_m"]

    #Take in the profile
    profile = block[profile_name,"rho_"+prof]
    
    #Calculate the 1 halo correlation function
    xi_1halo = calc_xi_1halo.calc_xi_1halo(profile,H0,OM)

    #Save to the block
    section = "xi_1halo"
    block.put_double_array_nd(section,"xi_1halo",xi_1halo)

    return 0
        
def cleanup(self):
    #No cleanup necessary
    return 0


