from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
import numpy as np

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
    return (r)

def execute(block,config):
    #Take in the configured parameters
    r = config

    #Write the masses and radii to the block
    section = "setup_delta_sigma"
    block[section,"radii"]=r
    return 0
        
def cleanup(self):
    #No cleanup necessary
    return 0
