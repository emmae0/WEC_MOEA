#### control file for running WAMIT with geometry defined by vector
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math
import os

### fixed parameters
pi = 3.14
g = 9.81

def general_WAMIT_control(geom_vector,make_mesh_figure,frequencies,NN,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55,beta_PTO,k_PTO,ulen,IALTFRC):

  [s1,r2,r1,h,l,a2,b2,s2] = geom_vector

  # geom_vector: s1,r2,r1,h,l,a2,b2,s2
  # make_mesh_figure: plot 2d and 3d plots of mesh
  # frequencies: which frequencies to run WAMIT for (list)
  # NN: number of panels on the smallest length
  # modes: surge sway heave roll pitch yaw (0 or 1)--which modes to run
  # betas: which angles for run WAMIT for (list)
  # water_depth: water depth
  # dir_name: directory path to run WAMIT (allows for parallel runs of WAMIT)
  # rho: water density
  # xCG,zCG: x and z coordinates of center of gravity
  # I55: pitch moment of inertia
  # beta_PTO and k_PTO: damping and stiffness for PTO
  
  # calculate volume of submerged body and mass of body (freely floating)
  from calculate_volume import calculate_volume
  Volume = calculate_volume(geom_vector)
  mass = Volume*rho
  
  ### delete old WAMIT files and change to correct directory
  from delete_WAMIT_files import delete_WAMIT_files
  delete_WAMIT_files(dir_name)
  
  ### make GDF, POT, FRC files
  from make_gdf import make_gdf
  from make_pot import make_pot
  if IALTFRC == 1:
    from make_frc_altform_1 import make_frc
    frc_string = make_frc(zCG,mass,I55,s2)
    
  elif IALTFRC == 2:
    from make_frc import make_frc
    frc_string = make_frc(rho,xCG,zCG,mass,I55,beta_PTO,k_PTO,s2)
  
  gdf_string = make_gdf(geom_vector,NN,make_mesh_figure,rho,ulen)
  pot_string = make_pot(water_depth,frequencies,betas,modes,s2)

  
  with open('mpswec.gdf','w+') as gdf_file:
    gdf_file.write(gdf_string)
  with open('mpswec.frc','w+') as frc_file:
    frc_file.write(frc_string)
  with open('mpswec.pot','w+') as pot_file:
    pot_file.write(pot_string)
  with open('fnames.wam','w+') as name_file:
    name_string = 'mpswec.pot \nmpswec.frc \nmpswec.cfg'
    name_file.write(name_string)
  with open('mpswec.cfg','w+') as cfg_file:
    cfg_string = '! MPS opt project \nILOWHI=0 \nIALTFRC =%s \nIPERIN=2 \nIPEROUT=2 \nIRR=0 \nILOG=1 \nISOR=1 \nMONITR=0 \nIDELFILES=4 \nISOLVE=0 '%str(IALTFRC)
    cfg_file.write(cfg_string)
    
  ### Run WAMIT
  os.system('wamit')
  
  #### open added mass and damping list
  
  A_B_out = list(open('mpswec.1','r+'))
  
    
  freq_list_out = []  
  A11_list = []
  A33_list = []
  A55_list = []
  B11_list = []
  B33_list = []
  B55_list = []
  
  
  for A_B_idx in range(len(A_B_out)):
    freq = float(A_B_out[A_B_idx][0:15])
    idx1 = float(A_B_out[A_B_idx][16:20])
    idx2 = float(A_B_out[A_B_idx][20:26])
    
    if A_B_idx == 0:
      freq_list_out.append(freq)
    else:  
      if freq_list_out[-1]!= freq:  
        freq_list_out.append(freq)
    if idx1 == 1 and idx2 == 1:
          A11_list.append(float(A_B_out[A_B_idx][27:41]))
          B11_list.append(float(A_B_out[A_B_idx][41:55])) 
    elif idx1 == 3 and idx2 == 3:
          A33_list.append(float(A_B_out[A_B_idx][27:41]))
          B33_list.append(float(A_B_out[A_B_idx][41:55]))
    elif idx1 == 5 and idx2 == 5:
          A55_list.append(float(A_B_out[A_B_idx][27:41]))
          B55_list.append(float(A_B_out[A_B_idx][41:55]))
      
  #freq_list_out = freq_list_out[::5]      
  
  ### Open .4 file (RAO's)
  xi_out = list(open('mpswec.4','r+'))
  
  xi_1_list = []
  xi_3_list = []
  xi_5_list = []
  
  for xi_idx in range(len(xi_out)):
    freq = float(xi_out[xi_idx][0:15])
    idx1 = float(xi_out[xi_idx][28:36])
    if idx1 == 1:
      xi_1_list.append(float(xi_out[xi_idx][36:49]))
    elif idx1 == 3:
      xi_3_list.append(float(xi_out[xi_idx][36:49]))
    elif idx1 == 5:
      #print (xi_out[xi_idx])
      #print (xi_out[xi_idx][36:49])
      #input()
      
      xi_5_list.append(float(xi_out[xi_idx][36:49]))
      
  ### Open .2 file (exciting force)
  X_out = list(open('mpswec.2','r+'))
  
  X_1_list = []
  X_3_list = []
  X_5_list = []
  
  for X_idx in range(len(X_out)):
    freq = float(X_out[X_idx][0:15])
    idx1 = float(X_out[X_idx][28:36])
    if idx1 == 1:
      X_1_list.append(float(X_out[X_idx][36:49]))
    elif idx1 == 3:
      X_3_list.append(float(X_out[X_idx][36:49]))
    elif idx1 == 5:
      X_5_list.append(float(X_out[X_idx][36:49]))
      
      
    
         
  return xi_1_list,xi_3_list,xi_5_list,A11_list,A33_list,A55_list,B11_list,B33_list,B55_list,freq_list_out,X_1_list,X_3_list,X_5_list,mass
    
  
  
#s1 = 3
#r2 = 3
#r1 = 1
#h = 4
#l = 6
#a2 = 0
#b2 = -0.5

#geom_vector_test = [s1,r2,r1,h,l,a2,b2]
#make_mesh_figure_test = 'yes'
#frequencies_test = [1]
#NN_test = 3
#modes_test = '1 0 1 0 1 0'
#betas_test = [0]
#water_depth_test = 50
#dir_name_test ='here'
#rho_test = 1000
#zCG_test = -1
#xCG_test = 1
#I55_test = 1
#beta_PTO_test = 5
#k_PTO_test = 0.5

#general_WAMIT_control(geom_vector_test,make_mesh_figure_test,frequencies_test,NN_test,modes_test,betas_test,water_depth_test,dir_name_test,rho_test,xCG_test,zCG_test,I55_test,beta_PTO_test,k_PTO_test)




