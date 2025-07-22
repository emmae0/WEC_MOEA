#### making FRC file for WAMIT for MPS WEC optimization study
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

def make_frc(zCG,mass,I55,zP):

  Header = 'WAMIT FRC file for MPS-UoP WEC optimization project'

  L_1 = '1 1 1 2 0 0 0 0 0' 
  # (Added mass and damping), (Exciting forces from Haskind), (Exciting forces from diffraction potential), (Motions of body), (Pressure/Velocity on body surface), (Pressure/free-surface elevation/Fluid velocity vector at field points), (Mean drift force from control surface), (Mean drift force from momentum), (Mean drift force from pressure)
  
  L_2 = '%s'%str(zCG-zP)
  
  ### mass matrix
  
  L_3 = '0 0 0' # set to 0 since no roll
  
  L_4 = '0 %s 0'%(str((I55/mass)**(1/2)))
  
  L_5 = '0 0 0' # set to 0 since no yaw
  
  ### angles
  
  L_6 = '0' # Haskind wave headings
  
  L_7 = '0' # Number of points where pressure and/or velocity are to be evaluated
  
  frc_string = "\n".join([Header,L_1,L_2,L_3,L_4,L_5,L_6,L_7])
  
  return frc_string
  
 
  
  
  
  
  
