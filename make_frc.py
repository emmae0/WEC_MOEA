#### making FRC file for WAMIT for MPS WEC optimization study
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

def make_frc(rho,xCG,zCG,mass,I55,beta_PTO,k_PTO,s2):


  Header = 'WAMIT FRC file for MPS-UoP WEC optimization project'

  L_1 = '1 1 1 -1 0 0 0 0 0' 
  # (Added mass and damping), (Exciting forces from Haskind), (Exciting forces from diffraction potential), (Motions of body), (Pressure/Velocity on body surface), (Pressure/free-surface elevation/Fluid velocity vector at field points), (Mean drift force from control surface), (Mean drift force from momentum), (Mean drift force from pressure)
  
  L_1_2 = '6'
  L_1_3 = '(0,0,0,0,1,0)'
  
  L_2 = str(rho)
  
  L_3 = '%s 0 %s'%(str(xCG),str(zCG-s2))
  #print (L_3)
  #input()
  
  ### mass matrix
 
  
  L_4 = '1' #yes, we want to specify mass matrix
  

  L_5 = '0 0 0 0 0 0'
  L_6 = '0 0 0 0 0 0'
  L_7 = '0 0 0 0 0 0'
  L_8 = '0 0 0 0 0 0'
  L_9 = '0 0 0 %s 0 0'%str(I55)
  L_10 = '0 0 0 0 0 0'
  
  
  #### damping matrix
  
  L_11 = '1' #yes, we want to specify external damping
  
  
  L_12 = '0 0 0 0 0 0 \n0 0 0 0 0 0 \n0 0 0 0 0 0 \n0 0 0 0 0 0 \n0 0 0 0 %s 0 \n0 0 0 0 0 0'%str(beta_PTO)
  
  #### stiffness matrix
  
  L_13 = '1' #yes, we want to specify external stiffness
  
  L_14 = '0 0 0 0 0 0 \n0 0 0 0 0 0 \n0 0 0 0 0 0 \n0 0 0 0 0 0 \n0 0 0 0 %s 0 \n0 0 0 0 0 0'%str(k_PTO)
  
  L_15 = '0' #no extra Haskind angles
  
  L_16 = '0' #no points where pressure/wave elevation/fluid velocity should be calculated
  

  frc_string = "\n".join([Header,L_1,L_1_2,L_1_3,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10,L_11,L_12,L_13,L_14,L_15,L_16])
  
  return frc_string
  
 
  
  
  
  
  
