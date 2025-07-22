#### making pot file for WAMIT for MPS WEC optimization
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

def make_pot(water_depth,frequencies,betas,modes,s2):

  Header = '! WAMIT pot file for MPS-UoP WEC optimization project'
  
  L_1 = str(water_depth)
  
  L_2 = '0 0' # this specifies how to solve the radiation and diffraction problem. I've specified that it only needs to solve the radiation and diffraction problems for the modes specified (not the 6 DOF problem)

  L_3 = str(len(frequencies))
  
  L_4 = ''
  for freq_idx in range(len(frequencies)):
    L_4 = L_4 + str(frequencies[freq_idx]) + ' '
    
  L_5 = str(len(betas))
  
  L_6 = ''
  for beta_idx in range(len(betas)):
    L_6 = L_6 + str(betas[beta_idx]) + ' '
    
  L_7 = '1' # number of bodies
  
  L_8 = 'mpswec.gdf' # name of gdf file
  
  L_9 = '0 0 %s 0'%str(s2) # coordinates and angle of body fixed coordinate system
  
  L_10 = modes
  
  pot_string = "\n".join([Header,L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10])
  
  return pot_string



