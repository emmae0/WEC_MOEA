### trying out multi-objective opt, varying b2 as test

import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import matplotlib.cm as cm

g = 9.81
pi = 3.14

### This is where we define the values that we want to consider in the optimisation.
# See paper for the definitions of these parameters
# Define them as lists
possible_l = 
possible_h = 
possible_s2 = 
possible_s1 = 

possible_r2 = 
possible_r1 = 
possible_a2 = 
possible_b2 = 

# Define filename directories for where to save the PF and pop. Use '.csv' at the end
PF_save_filename = '.csv'
pop_save_filename = '.csv'

# This is a flag for if you want to show plots to visualise the optimisation
show_opt_plots = 'no'

# vertical center of gravity.
zCG = 

# PTO spring (note: we always had it set to 0)
k_PTO = 0

### Set up WAMIT
# Which modes to allow (surge sway heave roll pitch yaw)
modes = '0 0 0 0 1 0'
# incident angles
betas = [0]
# water depth
water_depth = 
# directory for where the WAMIT files will be saved
dir_name = 
# water density
rho = 1026
# gravity
g = 9.81
# flag for if you want to visualise the mesh
make_mesh_figure = 'no'
# this is an integer showing the smallest number of nodes in a section for the first WAMIT wrapper function passthrough. Note that within the WAMIT wrapper, this is increased until it has converged (more notes in that function
NN_start = 5
# this is the percent error allowed to show convergence. That is: we require that adding one more panel to the smallest section results in a change of less than this percent
perc_error_wrapper = 0.04
# this is the absolute error allowed (note: it is an OR statement for convergence-- we say it converges if either the percent error is within the tolerance or the absolute error is within the tolerance)
abs_error_wrapper = 0.01
# this is the WAMIT length scale [X]
ulen = 1
# this is a WAMIT flag which allows us to specify that the body is only moving in pitch, centred about 
IALTFRC = 2

### incident wave period, frequency, wavenumber
T = 8
omega = 2*pi/T
k = 6.28189E-2

from MOEA import MOEA

### evolutionary algorithm parameters (populations size, number of generations, mutation probability)
popsize = 1500
num_gens = 4000
mut_prob = 0.05

### this runs the evolutionary algorithm - see that file for output definitions
pop_V_vectors,pop_F5,pop_one_over_kW,pop_geom_vectors,PF_V_vectors,PF_F5,PF_one_over_kW, PF_geom_vectors =  MOEA(make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k,popsize,mut_prob,num_gens,possible_s1,possible_r2,possible_r1,possible_h,possible_a2,possible_b2,possible_s2,possible_l,show_opt_plots)

print (PF_F5)
print ([1/w for w in PF_one_over_kW])
print (PF_geom_vectors)

# this plots the population and PF
plt.figure(0)
plt.scatter(pop_F5,[1/w for w in pop_one_over_kW],color='b')
plt.scatter(PF_F5,[1/w for w in PF_one_over_kW],color='g')
plt.xlabel('$F_5$')
plt.ylabel('$kW$')
plt.xlim(left=0)
plt.ylim(bottom=0)

plt.figure(1)
plt.scatter(pop_F5,pop_one_over_kW,color='b')
plt.scatter(PF_F5,PF_one_over_kW,color='g')
plt.xlabel('$F_5$')
plt.ylabel('$1/kW$')
plt.xlim(left=0)
plt.ylim(bottom=0)

from plot_2D import plot_2D

### this section plots the 2D shapes of the Pareto Front
colors = cm.rainbow(np.linspace(0,1,len(PF_F5)))

# order by increasing F5:
sorted_idx = np.argsort(np.asarray(PF_F5))
sorted_F5 = [PF_F5[idx] for idx in sorted_idx]
sorted_one_over_kW = [PF_one_over_kW[idx] for idx in sorted_idx]
sorted_geom_vectors = [PF_geom_vectors[idx] for idx in sorted_idx]

# this loops through the PF shapes and plots the 2D shapes
for PF_idx in range(len(PF_F5)):
  c = colors[PF_idx]
  geom_vector = sorted_geom_vectors[PF_idx]
  [s1,r2,r1,h,l,a2,b2,s2] = geom_vector
  
  x_c1,z_c1,x_c2,z_c2 = plot_2D(geom_vector)
  
  plt.figure(2)
  plt.plot(x_c1,z_c1,color=c)
  plt.plot(x_c2,z_c1,color=c)
  plt.plot([0,s1],[s2,0],color='k')
  plt.scatter([0],[s2],color='k',s=20)
  
  # this plots the Pareto Front with the colors corresponding to the 2D shapes
  plt.figure(3)
  plt.scatter(sorted_F5[PF_idx],[1/sorted_one_over_kW[PF_idx]],color=c)
  
### this saves the PF  
with open(PF_save_filename,'w') as csvfile:
  writer = csv.writer(csvfile,delimiter=' ')
  for PF_idx in range(len(PF_F5)):
    writer.writerow([PF_geom_vectors[PF_idx],PF_F5[PF_idx],1/PF_one_over_kW[PF_idx]])  
 
### this saves the population    
with open(pop_save_filename,'w') as csvfile:
  writer = csv.writer(csvfile,delimiter=' ')
  for pop_idx in range(len(pop_F5)):
    writer.writerow([pop_geom_vectors[pop_idx],pop_F5[pop_idx],1/pop_one_over_kW[pop_idx]])  

  
plt.figure(2)
plt.xlabel('x')
plt.ylabel('z')
plt.axhline(y=0,color='k',linestyle='--')
plt.axvline(x=0,color='k',linestyle='--')  

plt.figure(3)
plt.xlabel('$F_5$')
plt.ylabel('$kW$')
plt.xlim(left=0)
plt.ylim(bottom=0)
  

plt.show()

















