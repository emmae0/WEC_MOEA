### making text file of dimensions and mass properties of optimized shapes (Pareto Front)

import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import matplotlib.cm as cm
from matplotlib import gridspec

plt.rcParams["figure.autolayout"] = True 

g = 9.81
pi = 3.14


T_r = 8
omega_r = 2*pi/T_r
k_r = 6.28189E-2

omega_over_omega_r_list = np.concatenate((np.linspace(0.6,1.5,30)[0:12],np.linspace(0.95,1.05,10)[0:5],np.linspace(0.95,1.05,10)[5:],np.linspace(0.6,1.5,30)[15:]))
theta_list = np.linspace(0,90,10)


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



fig = plt.figure(figsize=(8,9))


gs = gridspec.GridSpec(3,1)#, width_ratios=[0.9],wspace=0.28, hspace=0.0, top=0.95, bottom=0.17, left=0.055, right=0.99) 
         
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[2,0])


from plot_2D import plot_2D

for kl_idx in range(3):

  if kl_idx == 0:
  
    PF_filename = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_0-5_IRR_1.csv'
    how_many_takeoff = 0
    marker_type = '*'
    markersize_type = 30
    l = 0.5/k_r
    linestyle_type = '--'
    #modes = '1 0 1 0 1 0'
  
  elif kl_idx == 1:

    PF_filename = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_1-1_IRR_1.csv'
    how_many_takeoff = 0
    marker_type = 'o'
    markersize_type = 45
    l = 18
    linestyle_type = '-'
    #modes = '0 0 0 0 1 0'
    
  elif kl_idx == 2:
  
    PF_filename = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_1-5_IRR_1.csv'
    how_many_takeoff = 0
    marker_type = '^'
    markersize_type = 30
    l = 1.5/k_r
    linestyle_type = 'dashdot'
    #modes = '1 0 1 0 1 0'
    
  PF_geom_vectors = []
  PF_one_over_kW = []
  PF_F5 = []
  with open(PF_filename) as csvfile:
    reader = csv.reader(csvfile,delimiter=',')
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
      
  colors = cm.rainbow(np.linspace(0,1,len(PF_F5)-how_many_takeoff))

  # order by increasing F5:
  sorted_idx = np.argsort(np.asarray(PF_F5))
  sorted_F5 = [PF_F5[idx] for idx in sorted_idx]
  sorted_one_over_kW = [PF_one_over_kW[idx] for idx in sorted_idx]
  sorted_geom_vectors = [PF_geom_vectors[idx] for idx in sorted_idx]
  

  for PF_idx in range(len(PF_F5)-how_many_takeoff):
    c = colors[PF_idx]
    geom_vector = sorted_geom_vectors[PF_idx]
    [s1,r2,r1,h,l,a2,b2,s2] = geom_vector
  
    x_c1,z_c1,x_c2,z_c2 = plot_2D(geom_vector)

    
    from calculate_xCB_and_zCB import calculate_xCB_and_zCB
    xCB,zCB = calculate_xCB_and_zCB(geom_vector)
   
    xCG = xCB
    
    ### first, calculate beta_55 = B_55 at omega_r and I_55
    
    from WAMIT_wrapper_func_IRR_1 import WAMIT_wrapper_func
    I55_in = 10
    beta_PTO_in = 0
    lookup = 'yes'
    A55_nd_r, B55_nd_r, X5_nd_r = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega_r],NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
  
    from calculate_volume import calculate_volume
    Volume = calculate_volume(geom_vector)
  
    S11 = l/3*((s1+r2)**3-(s1)**3)
    C55_1_nd = S11
    C55_2_nd = Volume*((zCB-s2)-(zCG-s2))
    C55 = rho*g*(C55_1_nd+C55_2_nd)
  
    A55_r = A55_nd_r*rho*ulen**5
    B55_r = B55_nd_r*rho*ulen**5*omega_r
    X5_r = X5_nd_r*rho*g*ulen**3
  
    I55 = (C55)/omega_r**2-A55_r

    mass = rho*Volume

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
    

    if kl_idx == 1:
      if PF_idx == 0:
        ax1.scatter([F5_nd_r],[xi_5_r*s1],color=c,marker=marker_type,s=markersize_type,edgecolors='k',label='$k_rl = %s$'%str(round(k_r*l,1)),zorder=100)
      else:
        ax1.scatter([F5_nd_r],[xi_5_r*s1],color=c,marker=marker_type,s=markersize_type,edgecolors='k',zorder=100)
    else:
      if PF_idx == 0:
        ax1.scatter([F5_nd_r],[xi_5_r*s1],color='grey',marker=marker_type,s=markersize_type,label='$k_rl = %s$'%str(round(k_r*l,1)))
      else:
        ax1.scatter([F5_nd_r],[xi_5_r*s1],color='grey',marker=marker_type,s=markersize_type)
      
    kW_list = []
    for omega_idx in range(len(omega_over_omega_r_list)):
    
      omega = round(omega_over_omega_r_list[omega_idx]*omega_r,3)
      k = omega**2/g

      lookup = 'yes'
        
      A55_nd, B55_nd, X5_nd = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega_r],NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
  
      A55 = A55_nd*rho*ulen**5
      B55 = B55_nd*rho*ulen**5*omega
      X5 = X5_nd*rho*g*ulen**3
  
      xi_5_sq = X5**2/((C55+k_PTO-omega**2*(I55+A55))**2+omega**2*(beta_PTO+B55)**2)
      xi_5 = xi_5_sq**(1/2)
  
      F5 = beta_PTO*omega*xi_5+abs(k_PTO)*xi_5
  
      P_5 = 0.5*beta_PTO*omega**2*xi_5**2
      
      Vg = omega/k*0.5*(1+2*k*water_depth/(math.sinh(2*k*water_depth)))

      P_I = 1/2*rho*g*Vg
  
      W_5 = P_5/P_I
  
      kW = k_r*W_5
      #kW = k*W_5
      F5_nd = F5/(rho*omega_r**2/k_r**4)

      kW_list.append(kW)
    
    
    ### half width at half height
    # find half max kW value
    half_max_kW = max(kW_list)/2
    max_kW_idx = kW_list.index(max(kW_list))
    # find closest index below
    first_half_kW = kW_list[:max_kW_idx]
    lower_idx = min(range(len(first_half_kW)), key=lambda i: abs(first_half_kW[i]-half_max_kW))
    #print (lower_idx)
    #print (half_max_kW)
    #print (kW_list[lower_idx])
    omega_lower = omega_over_omega_r_list[lower_idx]*omega_r
    second_half_kW = kW_list[max_kW_idx:]
    upper_idx = min(range(len(second_half_kW)), key=lambda i: abs(second_half_kW[i]-half_max_kW))+max_kW_idx
    #print (upper_idx)
    #print (kW_list[upper_idx])
    omega_upper = omega_over_omega_r_list[upper_idx]*omega_r
    full_width = omega_upper-omega_lower
    half_width = full_width/2
    
   
    if kl_idx == 1:
      if PF_idx == 0:
        ax2.scatter([F5_nd_r],[half_width],color=c,marker=marker_type,s=markersize_type,edgecolors='k',label='$k_rl = %s$'%str(round(k_r*l,1)),zorder=100)
      else:
        ax2.scatter([F5_nd_r],[half_width],color=c,marker=marker_type,s=markersize_type,edgecolors='k',zorder=100)
    else:
      if PF_idx == 0:
        ax2.scatter([F5_nd_r],[half_width],color='grey',marker=marker_type,s=markersize_type,label='$k_rl = %s$'%str(round(k_r*l,1)))
      else:
        ax2.scatter([F5_nd_r],[half_width],color='grey',marker=marker_type,s=markersize_type)





ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.set_xlabel('$|\widetilde{F_5}|$',fontsize=10)
ax1.set_ylabel('$|\\xi_5|s_1/A$',fontsize=10)
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.margins(x=0,y=0)
ax1.legend()


ax2.set_xlim(left=0)
ax2.set_ylim(bottom=0)
ax2.set_xlabel('$|\widetilde{F_5}|$',fontsize=10)
ax2.set_ylabel('$\Delta_{\omega} (rad)$',fontsize=10)
ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.margins(x=0,y=0)
ax2.legend()


modes = '0 0 0 0 1 0'

for kl_idx in range(3):

  if kl_idx == 0:
  
    PF_filename = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_0-5_IRR_1.csv'
    how_many_takeoff = 0
    marker_type = '+'
    markersize_type = 30
    l = 0.5/k_r
    linestyle_type = '--'
  
  elif kl_idx == 1:

    PF_filename = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_1-1_IRR_1.csv'
    how_many_takeoff = 0
    marker_type = 'o'
    markersize_type = 45
    l = 18
    linestyle_type = '-'
    
  elif kl_idx == 2:
  
    PF_filename = '/home/eedwards/MPS_WEC_opt/figures_and_tables_for_paper_IRR_1/7_r2_r1_b2_a2_kl_1-5_IRR_1.csv'
    how_many_takeoff = 0
    marker_type = '^'
    markersize_type = 30
    l = 1.5/k_r
    linestyle_type = 'dashdot'
    
  PF_geom_vectors = []
  PF_one_over_kW = []
  PF_F5 = []
  with open(PF_filename) as csvfile:
    reader = csv.reader(csvfile,delimiter=',')
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
      
  colors = cm.rainbow(np.linspace(0,1,len(PF_F5)-how_many_takeoff))

  # order by increasing F5:
  sorted_idx = np.argsort(np.asarray(PF_F5))
  sorted_F5 = [PF_F5[idx] for idx in sorted_idx]
  sorted_one_over_kW = [PF_one_over_kW[idx] for idx in sorted_idx]
  sorted_geom_vectors = [PF_geom_vectors[idx] for idx in sorted_idx]
  

  for PF_idx in range(len(PF_F5)-how_many_takeoff):

    c = colors[PF_idx]
    geom_vector = sorted_geom_vectors[PF_idx]
    [s1,r2,r1,h,l,a2,b2,s2] = geom_vector
  
    x_c1,z_c1,x_c2,z_c2 = plot_2D(geom_vector)
  

    from calculate_xCB_and_zCB import calculate_xCB_and_zCB
    xCB,zCB = calculate_xCB_and_zCB(geom_vector)
   
    xCG = xCB
    
    ### first, calculate beta_55 = B_55 at omega_r and I_55

    I55_in = 10
    beta_PTO_in = 0
    lookup = 'yes'
    A55_nd, B55_nd, X5_nd = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega_r],NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
    
    from calculate_volume import calculate_volume
    Volume = calculate_volume(geom_vector)
  
    S11 = l/3*((s1+r2)**3-(s1)**3)
    C55_1_nd = S11
    C55_2_nd = Volume*((zCB-s2)-(zCG-s2))
    C55 = rho*g*(C55_1_nd+C55_2_nd)
  
    A55_r = A55_nd*rho*ulen**5
    B55_r = B55_nd*rho*ulen**5*omega_r
    X5_r = X5_nd*rho*g*ulen**3
  
    I55 = (C55+k_PTO)/omega_r**2-A55_r

    mass = rho*Volume

    beta_PTO = B55_r
    
    xi_5_sq_r = X5_r**2/((C55+k_PTO-omega_r**2*(I55+A55_r))**2+omega_r**2*(beta_PTO+B55_r)**2)
    xi_5_r = xi_5_sq_r**(1/2)
  
    F5_r = beta_PTO*omega_r*xi_5_r+abs(k_PTO)*xi_5_r
    
    F5_nd_r = F5_r/(rho*omega_r**2*l/k_r**3)
    
    xi_5_list = []
    kW_list = []
    F5_list = []
    for theta_idx in range(len(theta_list)):
    
      theta = round(theta_list[theta_idx],1)
      
      lookup = 'yes'

      A55_nd, B55_nd, X5_nd = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega_r],NN_start,modes,[theta],water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
  
      A55 = A55_nd*rho*ulen**5
      B55 = B55_nd*rho*ulen**5*omega_r
      X5 = X5_nd*rho*g*ulen**3
  
      xi_5_sq = X5**2/((C55+k_PTO-omega_r**2*(I55+A55))**2+omega_r**2*(beta_PTO+B55)**2)
      xi_5 = xi_5_sq**(1/2)
  
      F5 = beta_PTO*omega_r*xi_5+abs(k_PTO)*xi_5
  
      P_5 = 0.5*beta_PTO*omega_r**2*xi_5**2
      
      Vg = omega_r/k_r*0.5*(1+2*k_r*water_depth/(math.sinh(2*k_r*water_depth)))

      P_I = 1/2*rho*g*Vg
  
      W_5 = P_5/P_I
  
      kW = k_r*W_5
      F5_nd = F5/(rho*omega_r**2/k_r**4)
      
      xi_5_list.append(xi_5)
      kW_list.append(kW)
      F5_list.append(F5_nd)
      
    
    
    ### half width at half height
    # find half max kW value
    half_max_kW = max(kW_list)/2
    max_kW_idx = kW_list.index(max(kW_list))
    # find closest index
    half_height_idx = min(range(len(kW_list)), key=lambda i: abs(kW_list[i]-half_max_kW))
    
    theta_half = theta_list[half_height_idx]*3.14/180
    
    
    
    if kl_idx == 1:
      if PF_idx == 0:
        ax3.scatter([F5_nd_r],[theta_half],color=c,marker=marker_type,s=markersize_type,edgecolors='k',label='$k_rl = %s$'%str(round(k_r*l,1)),zorder=100)
      else:
        ax3.scatter([F5_nd_r],[theta_half],color=c,marker=marker_type,s=markersize_type,edgecolors='k',zorder=100)
    else:
      if PF_idx == 0:
        ax3.scatter([F5_nd_r],[theta_half],color='grey',marker=marker_type,s=markersize_type,label='$k_rl = %s$'%str(round(k_r*l,1)))
      else:
        ax3.scatter([F5_nd_r],[theta_half],color='grey',marker=marker_type,s=markersize_type)
        
        

ax3.set_xlim(left=0)
ax3.set_ylim([0,1.3])
ax3.set_xlabel('$|\widetilde{F_5}|$',fontsize=10)
ax3.set_ylabel('$\Delta_{\\theta} (rad)$',fontsize=10)
ax3.tick_params(axis='both', which='major', labelsize=10)
ax3.legend()
ax3.margins(x=0,y=0)


#fig.subplots_adjust(wspace=0, hspace=0)

plt.show()
  
  

