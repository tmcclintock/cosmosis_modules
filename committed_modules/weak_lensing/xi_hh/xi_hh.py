from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
import numpy as np
import calc_xi_hh

#This module takes xi_mm and the bias^2 and multiplies them onto each other
def setup(options):
    #Take in the bias
    bias_choice = options[option_section,"bias"]
    return (bias_choice)

def execute(block,config):
    #Take in the configured parameters
    bias_choice = config+"_bias"#this is how the bias block is named

    #Acquire the bias
    bias = block[bias_choice,"bias"]

    #Acquire xi_mm
    xi_mm = block["xi_mm","xi_mm"]

    #Do the calculation
    xi_hh = calc_xi_hh.calc_xi_hh(bias,xi_mm)

    #Save to the block
    section = "xi_hh"
    block[section,"xi_hh"]=xi_hh

    return 0
        
def cleanup(self):
    #No cleanup necessary
    return 0


