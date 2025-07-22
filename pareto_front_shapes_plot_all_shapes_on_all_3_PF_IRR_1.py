### making text file of dimensions and mass properties of optimized shapes (Pareto Front)

import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import matplotlib.cm as cm

plt.rcParams["figure.autolayout"] = True    

g = 9.81
pi = 3.14
rho = 1026

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

PF_filename_1p1 = '/home/eedwards/MPS_WEC_opt/optimization_results/7_r2_r1_b2_a2_PF_corrected.csv'
PF_filename_0p5 = '/home/eedwards/MPS_WEC_opt/optimization_results/7_r2_r1_b2_a2_PF_corrected_kl_0-5.csv'
PF_filename_1p5 = '/home/eedwards/MPS_WEC_opt/optimization_results/7_r2_r1_b2_a2_PF_corrected_kl_1-5.csv'

save_PF_filename_1p1 = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_1-1_IRR_1.csv'
save_PF_filename_0p5 = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_0-5_IRR_1.csv'
save_PF_filename_1p5 = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_1-5_IRR_1.csv'

T_r = 8
omega_r = 2*pi/T_r
k_r = 6.28189E-2


zCG = -5

k_PTO = 0

modes = '0 0 0 0 1 0'
betas = [0]
water_depth = 85
dir_name = '/home/eedwards/MPS_WEC_opt/python_codes/WAMIT_files_1'
rho = 1026
g = 9.81
make_mesh_figure = 'no'
NN_start = 10
perc_error_wrapper = 0.04
abs_error_wrapper = 0.01
ulen = 1
IALTFRC = 2

l_1 = 18
l_2 = 0.5/k_r
l_3 = 1.5/k_r


PF_geom_vectors_1p1 = []
PF_one_over_kW_1p1 = []
PF_F5_1p1 = []
with open(PF_filename_1p1) as csvfile:
    reader = csv.reader(csvfile,delimiter=' ')
    for row in reader:
      geom_vector = row[0]
      gv_file_split = str(geom_vector).split(',')
      this_geom_vector_file = []
      for gv_idx in range(len(gv_file_split)):
        if gv_idx == 0:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][1:]))
        elif gv_idx == len(gv_file_split)-1:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][:-1]))
        else:
          this_geom_vector_file.append(float(gv_file_split[gv_idx]))
          
      PF_geom_vectors_1p1.append(this_geom_vector_file)
      PF_F5_1p1.append(float(row[1]))
      PF_one_over_kW_1p1.append(1/float(row[2]))
      
# order by increasing F5:
sorted_idx_1p1 = np.argsort(np.asarray(PF_F5_1p1))
sorted_F5_1p1 = [PF_F5_1p1[idx] for idx in sorted_idx_1p1]
sorted_one_over_kW_1p1 = [PF_one_over_kW_1p1[idx] for idx in sorted_idx_1p1]
sorted_geom_vectors_1p1 = [PF_geom_vectors_1p1[idx] for idx in sorted_idx_1p1] 

colors_1p1 = cm.rainbow(np.linspace(0,1,len(PF_F5_1p1)))

fig, ax1 = plt.subplots(figsize=(8,7))

for PF_idx in [len(PF_F5_1p1)-1]:#range(len(PF_F5_1p1)):
    c = colors_1p1[PF_idx]

    
    geom_vector = sorted_geom_vectors_1p1[PF_idx]
    [s1,r2,r1,h,l,a2,b2,s2] = geom_vector
  
    from calculate_xCB_and_zCB import calculate_xCB_and_zCB
    xCB,zCB = calculate_xCB_and_zCB(geom_vector)
    
    print (xCB, zCB)
    input()
   
    xCG = xCB
    
    from WAMIT_wrapper_func_IRR_1 import WAMIT_wrapper_func
    I55_in = 10
    beta_PTO_in = 0
    lookup = 'yes'
    A55_nd_IRR_1, B55_nd_IRR_1, X5_nd_IRR_1 = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega_r],NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
  
    from calculate_volume import calculate_volume
    Volume = calculate_volume(geom_vector)  
  
    S11 = l/3*((s1+r2)**3-(s1)**3)
    C55_1_nd = S11
    C55_2_nd = Volume*((zCB-s2)-(zCG-s2))
    C55 = rho*g*(C55_1_nd+C55_2_nd)
  
    A55_r = A55_nd_IRR_1*rho*ulen**5
    B55_r = B55_nd_IRR_1*rho*ulen**5*omega_r
    X5_r = X5_nd_IRR_1*rho*g*ulen**3
  
    I55 = (C55)/omega_r**2-A55_r

    beta_PTO = B55_r
    
    print (A55_r)
    print (B55_r)
    print (X5_r)
    print (C55)
    print (I55)
    print (C55+k_PTO-omega_r**2*(I55+A55_r))
    input()
    
    xi_5_sq_r = X5_r**2/((C55+k_PTO-omega_r**2*(I55+A55_r))**2+omega_r**2*(beta_PTO+B55_r)**2)
    xi_5_r = xi_5_sq_r**(1/2)
  
    F5_r = beta_PTO*omega_r*xi_5_r+abs(k_PTO)*xi_5_r
  
    P_5_r = 0.5*beta_PTO*omega_r**2*xi_5_r**2
      
    Vg_r = 0.5*omega_r/k_r#*(1+2*k*water_depth/(math.sinh(2*k*water_depth)))

    P_I_r = 1/2*rho*g*Vg_r
  
    W_5_r = P_5_r/P_I_r
  
    kW_r = k_r*W_5_r
    F5_nd_r = F5_r/(rho*omega_r**2*l/k_r**3)
    
    print (kW_r, F5_nd_r)
    input()
  
    plt.figure(0)
    plt.scatter([F5_nd_r],[kW_r],color=c,edgecolors='k',s=45)
    
    #F5_new_nd = sorted_F5_1p1[PF_idx]*rho*omega_r**2/k_r**4/(rho*omega_r**2*l_1/k_r**3)
    
    if PF_idx == 0:
      ax1.scatter([F5_nd_r],[kW_r],color=c,edgecolors='k',s=45,label='$k_r l = 1.1$',zorder=100)
    else:
      ax1.scatter([F5_nd_r],[kW_r],color=c,edgecolors='k',s=45,zorder=100)
    
    with open(save_PF_filename_1p1, 'a') as f:
      writer = csv.writer(f)
      writer.writerow([geom_vector,F5_nd_r,kW_r])
    
      
PF_geom_vectors_0p5 = []
PF_one_over_kW_0p5 = []
PF_F5_0p5 = []
with open(PF_filename_0p5) as csvfile:
    reader = csv.reader(csvfile,delimiter=' ')
    for row in reader:
      geom_vector = row[0]
      gv_file_split = str(geom_vector).split(',')
      this_geom_vector_file = []
      for gv_idx in range(len(gv_file_split)):
        if gv_idx == 0:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][1:]))
        elif gv_idx == len(gv_file_split)-1:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][:-1]))
        else:
          this_geom_vector_file.append(float(gv_file_split[gv_idx]))
          
      PF_geom_vectors_0p5.append(this_geom_vector_file)
      PF_F5_0p5.append(float(row[1]))
      PF_one_over_kW_0p5.append(1/float(row[2]))


# order by increasing F5:
sorted_idx_0p5 = np.argsort(np.asarray(PF_F5_0p5))
sorted_F5_0p5 = [PF_F5_0p5[idx] for idx in sorted_idx_0p5]
sorted_one_over_kW_0p5 = [PF_one_over_kW_0p5[idx] for idx in sorted_idx_0p5]
sorted_geom_vectors_0p5 = [PF_geom_vectors_0p5[idx] for idx in sorted_idx_0p5]

colors_0p5 = cm.rainbow(np.linspace(0,1,len(PF_F5_0p5)))

for PF_idx in range(len(PF_F5_0p5)):
    c = colors_0p5[PF_idx]
    
    geom_vector = sorted_geom_vectors_0p5[PF_idx]
    [s1,r2,r1,h,l,a2,b2,s2] = geom_vector
  
    from calculate_xCB_and_zCB import calculate_xCB_and_zCB
    xCB,zCB = calculate_xCB_and_zCB(geom_vector)
   
    xCG = xCB
    
    from WAMIT_wrapper_func_IRR_1 import WAMIT_wrapper_func
    I55_in = 10
    beta_PTO_in = 0
    lookup = 'yes'
    A55_nd_IRR_1, B55_nd_IRR_1, X5_nd_IRR_1 = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega_r],NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
  
    from calculate_volume import calculate_volume
    Volume = calculate_volume(geom_vector)
  
    S11 = l/3*((s1+r2)**3-(s1)**3)
    C55_1_nd = S11
    C55_2_nd = Volume*((zCB-s2)-(zCG-s2))
    C55 = rho*g*(C55_1_nd+C55_2_nd)
  
    A55_r = A55_nd_IRR_1*rho*ulen**5
    B55_r = B55_nd_IRR_1*rho*ulen**5*omega_r
    X5_r = X5_nd_IRR_1*rho*g*ulen**3
  
    I55 = (C55)/omega_r**2-A55_r

    beta_PTO = B55_r
    
    xi_5_sq_r = X5_r**2/((C55+k_PTO-omega_r**2*(I55+A55_r))**2+omega_r**2*(beta_PTO+B55_r)**2)
    xi_5_r = xi_5_sq_r**(1/2)
  
    F5_r = beta_PTO*omega_r*xi_5_r+abs(k_PTO)*xi_5_r
  
    P_5_r = 0.5*beta_PTO*omega_r**2*xi_5_r**2
      
    Vg_r = 0.5*omega_r/k_r#*(1+2*k*water_depth/(math.sinh(2*k*water_depth)))

    P_I_r = 1/2*rho*g*Vg_r
  
    W_5_r = P_5_r/P_I_r
  
    kW_r = k_r*W_5_r
    F5_nd_r = F5_r/(rho*omega_r**2*l/k_r**3)
      
    plt.figure(0)
    plt.scatter([F5_nd_r],[kW_r],color='grey',marker='+',s=30)
    
    if PF_idx == 0:
      ax1.scatter([F5_nd_r],[kW_r],color='grey',marker='+',s=30,label='$k_rl = 0.5$')
    else:
      ax1.scatter([F5_nd_r],[kW_r],color='grey',marker='+',s=30)
    
    with open(save_PF_filename_0p5, 'a') as f:
      writer = csv.writer(f)
      writer.writerow([geom_vector,F5_nd_r,kW_r])
    
    
PF_geom_vectors_1p5 = []
PF_one_over_kW_1p5 = []
PF_F5_1p5 = []
with open(PF_filename_1p5) as csvfile:
    reader = csv.reader(csvfile,delimiter=' ')
    for row in reader:
      geom_vector = row[0]
      gv_file_split = str(geom_vector).split(',')
      this_geom_vector_file = []
      for gv_idx in range(len(gv_file_split)):
        if gv_idx == 0:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][1:]))
        elif gv_idx == len(gv_file_split)-1:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][:-1]))
        else:
          this_geom_vector_file.append(float(gv_file_split[gv_idx]))
          
      PF_geom_vectors_1p5.append(this_geom_vector_file)
      PF_F5_1p5.append(float(row[1]))
      PF_one_over_kW_1p5.append(1/float(row[2]))


# order by increasing F5:
sorted_idx_1p5 = np.argsort(np.asarray(PF_F5_1p5))
sorted_F5_1p5 = [PF_F5_1p5[idx] for idx in sorted_idx_1p5]
sorted_one_over_kW_1p5 = [PF_one_over_kW_1p5[idx] for idx in sorted_idx_1p5]
sorted_geom_vectors_1p5 = [PF_geom_vectors_1p5[idx] for idx in sorted_idx_1p5]

colors_1p5 = cm.rainbow(np.linspace(0,1,len(PF_F5_1p5)-5))

for PF_idx in range(len(PF_F5_1p5)-5):
    c = colors_1p5[PF_idx]
    
    geom_vector = sorted_geom_vectors_1p5[PF_idx]
    [s1,r2,r1,h,l,a2,b2,s2] = geom_vector
  
    from calculate_xCB_and_zCB import calculate_xCB_and_zCB
    xCB,zCB = calculate_xCB_and_zCB(geom_vector)
   
    xCG = xCB
    
    from WAMIT_wrapper_func_IRR_1 import WAMIT_wrapper_func
    I55_in = 10
    beta_PTO_in = 0
    lookup = 'yes'
    A55_nd_IRR_1, B55_nd_IRR_1, X5_nd_IRR_1 = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega_r],NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
  
    from calculate_volume import calculate_volume
    Volume = calculate_volume(geom_vector)
  
    S11 = l/3*((s1+r2)**3-(s1)**3)
    C55_1_nd = S11
    C55_2_nd = Volume*((zCB-s2)-(zCG-s2))
    C55 = rho*g*(C55_1_nd+C55_2_nd)
  
    A55_r = A55_nd_IRR_1*rho*ulen**5
    B55_r = B55_nd_IRR_1*rho*ulen**5*omega_r
    X5_r = X5_nd_IRR_1*rho*g*ulen**3
  
    I55 = (C55)/omega_r**2-A55_r

    beta_PTO = B55_r
    
    xi_5_sq_r = X5_r**2/((C55+k_PTO-omega_r**2*(I55+A55_r))**2+omega_r**2*(beta_PTO+B55_r)**2)
    xi_5_r = xi_5_sq_r**(1/2)
  
    F5_r = beta_PTO*omega_r*xi_5_r+abs(k_PTO)*xi_5_r
  
    P_5_r = 0.5*beta_PTO*omega_r**2*xi_5_r**2
      
    Vg_r = 0.5*omega_r/k_r#*(1+2*k*water_depth/(math.sinh(2*k*water_depth)))

    P_I_r = 1/2*rho*g*Vg_r
  
    W_5_r = P_5_r/P_I_r
  
    kW_r = k_r*W_5_r
    F5_nd_r = F5_r/(rho*omega_r**2*l/k_r**3)
  
    plt.figure(0)
    plt.scatter([F5_nd_r],[kW_r],color='r',marker='^',s=30)
        
    
    if PF_idx == 0:
      ax1.scatter([F5_nd_r],[kW_r],color='grey',marker='^',s=30,label='$k_rl = 1.5$')
    else:
      ax1.scatter([F5_nd_r],[kW_r],color='grey',marker='^',s=30)
    ax1.legend()

    with open(save_PF_filename_1p5, 'a') as f:
      writer = csv.writer(f)
      writer.writerow([geom_vector,F5_nd_r,kW_r])
    
    

plt.figure(0)
plt.xlim(left=0)
plt.ylim(bottom=0)      
plt.xlabel('$F_5 k_r^4/(\\rho \omega_r^2 A)$',fontsize=10)
plt.ylabel('$k_r W$',fontsize=10)

plt.figure(1)
plt.xlim(left=0)
plt.ylim(bottom=0)      
plt.xlabel('$|\widetilde{F_5}|$',fontsize=10)
plt.ylabel('$k_r W$',fontsize=10)


PF_filename = '/home/eedwards/MPS_WEC_opt/optimization_results/7_r2_r1_b2_a2_PF_corrected.csv'
how_many_takeoff = 0
l = 18

T_r = 8
omega_r = 2*pi/T_r
k_r = 6.28189E-2


PF_geom_vectors = []
PF_one_over_kW = []
PF_F5 = []
with open(PF_filename) as csvfile:
    reader = csv.reader(csvfile,delimiter=' ')
    for row in reader:
      geom_vector = row[0]
      gv_file_split = str(geom_vector).split(',')
      this_geom_vector_file = []
      for gv_idx in range(len(gv_file_split)):
        if gv_idx == 0:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][1:]))
        elif gv_idx == len(gv_file_split)-1:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][:-1]))
        else:
          this_geom_vector_file.append(float(gv_file_split[gv_idx]))
          
      PF_geom_vectors.append(this_geom_vector_file)
      PF_F5.append(float(row[1]))
      PF_one_over_kW.append(1/float(row[2]))


zCG = -5

k_PTO = 0

modes = '0 0 0 0 1 0'
betas = [0]
water_depth = 85
dir_name = '/home/eedwards/MPS_WEC_opt/python_codes/WAMIT_files_3'
rho = 1026
g = 9.81
make_mesh_figure = 'no'
NN_start = 19
perc_error_wrapper = 0.04
abs_error_wrapper = 0.01
ulen = 1
IALTFRC = 2




colors = cm.rainbow(np.linspace(0,1,len(PF_F5)-how_many_takeoff))

# order by increasing F5:
sorted_idx = np.argsort(np.asarray(PF_F5))
sorted_F5 = [PF_F5[idx] for idx in sorted_idx]
sorted_one_over_kW = [PF_one_over_kW[idx] for idx in sorted_idx]
sorted_geom_vectors = [PF_geom_vectors[idx] for idx in sorted_idx]

F5_normalized_list = [f/sorted_F5[-1] for f in sorted_F5]
kW_normalized_list = [(1/w)/(1/sorted_one_over_kW[-1]) for w in sorted_one_over_kW]


from plot_2D import plot_2D


for PF_idx in range(len(PF_F5)-how_many_takeoff):

    c = colors[PF_idx]

    F5_dim = sorted_F5[PF_idx]*rho*omega_r**2/k_r**4
    F5_nd = F5_dim/(rho*omega_r**2*l/k_r**3)

    geom_vector = sorted_geom_vectors[PF_idx]
    [s1,r2,r1,h,l,a2,b2,s2] = geom_vector
    
    
    x_c1,z_c1,x_c2,z_c2 = plot_2D(geom_vector)

    kx_c1 = [round(k_r*xx,3) for xx in x_c1]
    kx_c2 = [round(k_r*xx,3) for xx in x_c2]
    kz_c1 = [round(k_r*xx,3) for xx in z_c1]
    kz_c2 = [round(k_r*xx,3) for xx in z_c2]
    
    
    
    width, height = [0.13, 0.13] #(1.4,1.2)
    left_offset = 0.04
    
    if PF_idx < 5:
      left = (PF_idx+1)*width+left_offset
      bottom = 0.41
      ax2 = fig.add_axes([left, bottom, width, height])
      ax2.tick_params(axis='x',which='both',labelsize=8) 
    elif PF_idx < 11:
      left = (PF_idx-4)*width+left_offset
      bottom = 0.28
      ax2 = fig.add_axes([left, bottom, width, height])
      ax2.tick_params(axis='x',which='both',labelsize=8) 
    else:
      left = (PF_idx-10)*width+left_offset
      bottom = 0.15
      ax2 = fig.add_axes([left, bottom, width, height])
      ax2.tick_params(axis='x',which='both',labelsize=8) 
    
    
    ax2.tick_params(axis='y',which='both',labelleft=False,labelsize=8)  
     

    ax2.plot(kx_c1,kz_c1,color='k',linewidth=1)
    ax2.plot(kx_c2,kz_c2,color='k',linewidth=1)
    if a2 >= 0 or r1 >= 0:
      ax2.fill_between(kx_c1,kz_c1,color=c)
    ax2.fill_between(kx_c2,kz_c2,color=c)
    
    in_both_kx = []
    in_both_kz1 = []
    in_both_kz2 = []
    for c_idx in range(len(kx_c1)):
      if kx_c1[c_idx]>=kx_c2[-1]+0.001:
          in_both_kx.append(kx_c1[c_idx])
          in_both_kz1.append(kz_c1[c_idx])

          in_both_kz2.append(kz_c2[find_nearest(kx_c2,kx_c1[c_idx])])
        
    
   
    
    if r1 < s1:
      ax2.fill_between(in_both_kx,in_both_kz1,color='white')
      
    else:
      ax2.fill_between(in_both_kx,in_both_kz2,in_both_kz1,color='white')
      
    print (s1,s1+r2)

    ax2.axhline(y=0,color='white')
    ax2.plot([0,k_r*s1],[k_r*s2,0],color='k',linewidth=1)
    ax2.scatter([0],[k_r*s2],color='k',s=20)
    ax2.plot([k_r*s1,k_r*(s1+r2)],[0,0],color='k',linewidth=1)
    if PF_idx > 10:
      ax2.set_xlabel('$k_rx$',fontsize=8)
    if PF_idx == 0 or PF_idx == 5 or PF_idx == 11:
      ax2.set_ylabel('$k_rz$',fontsize=8)
      ax2.tick_params(axis='both', which='major', labelsize=8,labelleft=True)
    ax2.axhline(y=0,color='k',linestyle='--',linewidth=1)
    ax2.axvline(x=0,color='k',linestyle='--',linewidth=1)  
    ax2.set_xlim([-2*k_r,k_r*15])

    
  
    #plt.show()
   

    

plt.figure(0)
plt.xlim(left=0)
plt.ylim(bottom=0)      
plt.xlabel('$|\widetilde{F_5}|$',fontsize=10)
plt.ylabel('$k_rW$',fontsize=10)


  
plt.show()
  
  
  

