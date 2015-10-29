from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
import numpy as np

cosmo = section_names.cosmological_parameters
distances = section_names.distances
def setup(options):
    return (0)

#This module combines the 1 halo and 2 halo correlation functions
#into the cluster-matter corrlation function as per Hayashi and White
def execute(block,config):
    #Acquire the 1 and 2 halo correlation functions
    xi_1halo = block["xi_1halo","xi_1halo"]
    xi_2halo = block["xi_2halo","xi_2halo"]

    #Take the maximum of the two arrays
    #Also calculate the reduced arrays
    mask1 = xi_1halo>xi_2halo
    mask2 = xi_1halo<xi_2halo
    xi = xi_1halo*mask1

    xi+= xi_2halo*mask2
    xi_1halo_with_exclusion = xi_1halo*mask1
    xi_1halo_with_exclusion+= 0*mask2
    xi_2halo_with_exclusion = xi_2halo*mask2
    xi_2halo_with_exclusion+= 0*mask1

    section = "xi_hm"
    
    #Write all of the xi to the block
    block[section,"xi_hm"]=xi
    block[section,"xi_1halo_no_exclusion"]=xi_1halo
    block[section,"xi_2halo_no_exclusion"]=xi_2halo
    block[section,"xi_1halo_with_exclusion"]=xi_1halo_with_exclusion
    block[section,"xi_2halo_with_exclusion"]=xi_2halo_with_exclusion

    return 0
        
def cleanup(self):
    #No cleanup necessary
    return 0


