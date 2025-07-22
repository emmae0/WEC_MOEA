### trying out multi-objective opt, varying b2 as test

import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import matplotlib.cm as cm

g = 9.81
pi = 3.14

T = 8
omega = 2*pi/T
k = 6.28189E-2

kh = 0.2


possible_l = [18]
possible_h = [round(kh/k,2)]
possible_s2 = [4]
possible_s1 = [2.587]

possible_r2 = np.linspace(2,10,7)
possible_r1 = np.linspace(-7.4,12.587,7)
possible_a2 = np.linspace(-1,1,7)
possible_b2 = np.linspace(-4,1,7)


PF_save_filename = '/home/eedwards/MPS_WEC_opt/optimization_results/7_r2_r1_b2_a2_PF_corrected_kh_0-2.csv'
pop_save_filename = '/home/eedwards/MPS_WEC_opt/optimization_results/7_r2_r1_b2_a2_pop_corrected_kh_0-2.csv'

#[s1,r2,r1,h,l,a2,b2,s2]

show_opt_plots = 'no'

zCG = -5

k_PTO = 0

modes = '0 0 0 0 1 0'
betas = [0]
water_depth = 85
dir_name = '/home/eedwards/MPS_WEC_opt/python_codes/WAMIT_files_3'
rho = 1026
g = 9.81
make_mesh_figure = 'no'
NN_start = 5
perc_error_wrapper = 0.04
abs_error_wrapper = 0.01
ulen = 1
IALTFRC = 2



from MOEA import MOEA

popsize = 360
num_gens = 960
mut_prob = 0.05


pop_V_vectors,pop_F5,pop_one_over_kW,pop_geom_vectors,PF_V_vectors,PF_F5,PF_one_over_kW, PF_geom_vectors =  MOEA(make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k,popsize,mut_prob,num_gens,possible_s1,possible_r2,possible_r1,possible_h,possible_a2,possible_b2,possible_s2,possible_l,show_opt_plots)

print (PF_F5)
print ([1/w for w in PF_one_over_kW])
print (PF_geom_vectors)


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


colors = cm.rainbow(np.linspace(0,1,len(PF_F5)))

# order by increasing F5:
sorted_idx = np.argsort(np.asarray(PF_F5))
sorted_F5 = [PF_F5[idx] for idx in sorted_idx]
sorted_one_over_kW = [PF_one_over_kW[idx] for idx in sorted_idx]
sorted_geom_vectors = [PF_geom_vectors[idx] for idx in sorted_idx]

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
  
  plt.figure(3)
  plt.scatter(sorted_F5[PF_idx],[1/sorted_one_over_kW[PF_idx]],color=c)
  
with open(PF_save_filename,'w') as csvfile:
  writer = csv.writer(csvfile,delimiter=' ')
  for PF_idx in range(len(PF_F5)):
    writer.writerow([PF_geom_vectors[PF_idx],PF_F5[PF_idx],1/PF_one_over_kW[PF_idx]])  
    
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

















