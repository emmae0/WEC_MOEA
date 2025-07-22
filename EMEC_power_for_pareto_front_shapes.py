### making text file of dimensions and mass properties of optimized shapes (Pareto Front)

import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import matplotlib.cm as cm
from matplotlib import gridspec
import scipy

plt.rcParams["figure.autolayout"] = True 

g = 9.81
pi = 3.14

def spectrum(peak_per,swh,per):
  omega_m = 2*pi/peak_per
  omega = 2*pi/per
  f = 1/per
  A_s = 173*swh**2/peak_per**4
  B_s = 691/peak_per**4
  return A_s/((2*pi)**4*f**5)*math.exp(-B_s/((2*pi)**4*f**4))

#### occurence matrix
EMEC_occ_period = np.arange(0.5,19.5,1)
EMEC_occ_swh = np.arange(0.25,11.75,0.5)
EMEC_occ_matrix = np.zeros((len(EMEC_occ_swh),len(EMEC_occ_period)))
power_matrix = np.zeros((len(EMEC_occ_swh),len(EMEC_occ_period)))
force_matrix = np.zeros((len(EMEC_occ_swh),len(EMEC_occ_period)))
power_times_occ_matrix = np.zeros((len(EMEC_occ_swh),len(EMEC_occ_period)))
force_times_occ_matrix = np.zeros((len(EMEC_occ_swh),len(EMEC_occ_period)))
P_over_J_matrix = np.zeros((len(EMEC_occ_swh),len(EMEC_occ_period)))
k_times_P_over_J_matrix = np.zeros((len(EMEC_occ_swh),len(EMEC_occ_period)))
occ_times_k_times_P_over_J_matrix = np.zeros((len(EMEC_occ_swh),len(EMEC_occ_period)))


EMEC_occ_matrix[0,5:7] = np.matrix([0.1,0.1])/100
EMEC_occ_matrix[1,2:13] = np.matrix([0.01,0.04,0.16,0.28,0.33,0.23,0.13,0.07,0.03,0.01,0.01])/100
EMEC_occ_matrix[2,3:15] = np.matrix([0.05,0.3,0.68,0.95,0.94,0.7,0.41,0.2,0.09,0.04,0.02,0.01])/100
EMEC_occ_matrix[3,3:17] = np.matrix([0.01,0.22,0.83,1.25,1.49,1.25,0.88,0.52,0.26,0.11,0.06,0.02,0.01,0.01])/100
EMEC_occ_matrix[4,4:17] = np.matrix([0.07,0.79,1.71,2.33,2.19,1.36,0.8,0.42,0.19,0.08,0.04,0.03,0.02])/100
EMEC_occ_matrix[5,4:16] = np.matrix([0.01,0.43,1.7,2.64,2.5,1.72,0.98,0.52,0.26,0.12,0.04,0.02])/100
EMEC_occ_matrix[6,5:16] = np.matrix([0.15,1.31,2.71,2.97,2.19,1.34,0.65,0.31,0.17,0.09,0.03])/100
EMEC_occ_matrix[7,5:16] = np.matrix([0.05,0.78,2.65,2.94,2.12,1.37,0.64,0.31,0.11,0.04,0.04])/100
EMEC_occ_matrix[8,6:15] = np.matrix([0.33,1.91,2.89,2.31,1.39,0.7,0.27,0.1,0.05])/100
EMEC_occ_matrix[9,6:14] = np.matrix([0.09,1.17,2.73,2.46,1.41,0.72,0.22,0.06])/100
EMEC_occ_matrix[10,7:14] = np.matrix([0.7,2,2.13,1.38,0.57,0.14,0.07])/100
EMEC_occ_matrix[11,7:13] = np.matrix([0.2,1.34,2.24,1.1,0.45,0.16])/100
EMEC_occ_matrix[12,7:13] = np.matrix([0.06,0.72,1.69,1.3,0.36,0.1])/100
EMEC_occ_matrix[13,8:13] = np.matrix([0.31,0.94,0.85,0.31,0.11])/100
EMEC_occ_matrix[14,8:13] = np.matrix([0.18,0.49,0.66,0.36,0.13])/100
EMEC_occ_matrix[15,8:13] = np.matrix([0.1,0.34,0.62,0.27,0.15])/100
EMEC_occ_matrix[16,9:13] = np.matrix([0.13,0.28,0.31,0.17])/100
EMEC_occ_matrix[17,9:13] = np.matrix([0.14,0.16,0.35,0.19])/100
EMEC_occ_matrix[18,10:12] = np.matrix([0.18,0.19])/100


T_r = 8
omega_r = 2*pi/T_r
k_r = 6.28189E-2

omega_over_omega_r_list = np.concatenate((np.linspace(0.6,1.5,30)[0:12],np.linspace(0.95,1.05,10)[0:5],np.linspace(0.95,1.05,10)[5:],np.linspace(0.6,1.5,30)[15:]))

omega_list = [o*omega_r for o in omega_over_omega_r_list]
T_list = [2*pi/o for o in omega_list]

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


from plot_2D import plot_2D

PF_filename = '/home/eedwards/MPS_WEC_opt/optimization_results/7_r2_r1_b2_a2_PF_corrected.csv'
l = 18
    
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
      
# order by increasing F5:
sorted_idx = np.argsort(np.asarray(PF_F5))
sorted_F5 = [PF_F5[idx] for idx in sorted_idx]
sorted_one_over_kW = [PF_one_over_kW[idx] for idx in sorted_idx]
sorted_geom_vectors = [PF_geom_vectors[idx] for idx in sorted_idx]


power_list = []
max_power_list = []
force_list = []
max_force_list = []
k_P_over_J_list = []

for PF_idx in [8]:#range(len(sorted_F5)):
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
    
    for per_idx in [10]:#range(len(EMEC_occ_period)):
      for A_idx in [10]:#range(len(EMEC_occ_swh)):
        T_peak = EMEC_occ_period[per_idx]
        omega_peak = 2*pi/T_peak
        
        SWH = (EMEC_occ_swh[A_idx])
        
        P_list = []
        S_list = []
        P_times_S_list = [] #times 2
        f_list = []
        F_list = []
        F_times_S_list = []
        for omega in omega_list:
          k = omega**2/g

          lookup = 'yes'
          
          f = omega/(2*pi)
          f_list.append(f)
        
          A55_nd, B55_nd, X5_nd = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega_r],NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
  
          A55 = A55_nd*rho*ulen**5
          B55 = B55_nd*rho*ulen**5*omega
          X5 = X5_nd*rho*g*ulen**3
  
          xi_5_sq = X5**2/((C55+k_PTO-omega**2*(I55+A55))**2+omega**2*(beta_PTO+B55)**2)
          xi_5 = xi_5_sq**(1/2)
  
          F5 = beta_PTO*omega*xi_5+abs(k_PTO)*xi_5
  
          P_5 = 0.5*beta_PTO*omega**2*xi_5**2
          
          print (T_peak,SWH,2*pi/omega)
          S = spectrum(T_peak,SWH,2*pi/omega)
          
          P_list.append(P_5)
          S_list.append(S)
          P_times_S_list.append(2*P_5*S)
          
          F_list.append(F5)
          F_times_S_list.append((2*S)**(1/2)*F5)
          
        total_power = scipy.integrate.trapz(P_times_S_list,f_list)     
        print ('total power:',total_power)
        print ('total power in MW:',total_power/1e6)
        
        total_force = scipy.integrate.trapz(F_times_S_list,f_list)
        
        plt.figure(0)
        plt.plot(f_list,S_list)
        plt.xlabel('f')
        plt.ylabel('S')
        
        plt.figure(1)
        plt.plot(f_list,P_times_S_list)
        plt.xlabel('f')
        plt.ylabel('P*S')
        
        print (scipy.integrate.trapz(S_list,f_list))
        plt.show()
        
        #plt.figure(0)
        #plt.plot(f_list,S_list)
        
        #plt.figure(1)
        #plt.plot(f_list,P_times_S_list)
        
        #plt.figure(2)
        #plt.plot(f_list,[p/1e6 for p in P_list])
            
          
          
        power_matrix[A_idx,per_idx] = total_power
        power_times_occ_matrix[A_idx,per_idx] = EMEC_occ_matrix[A_idx,per_idx]*total_power*(EMEC_occ_period[1]-EMEC_occ_period[0])
        force_matrix[A_idx,per_idx] = total_force
        force_times_occ_matrix[A_idx,per_idx] = EMEC_occ_matrix[A_idx,per_idx]*total_force*(EMEC_occ_period[1]-EMEC_occ_period[0])

          
    print ('shape:',PF_idx)
    print ('total power in MW:',np.sum(power_times_occ_matrix)/1e6)
    print ('max power in MW:',np.max(power_matrix)/1e6)
    print ('max force in MN:',np.max(force_matrix)/1e6)
    print ('mean force in MN:',np.sum(force_times_occ_matrix)/1e6)
    print (np.argmax(force_matrix))
    force_max_idx = np.unravel_index(np.argmax(force_matrix, axis=None), force_matrix.shape)
    force_max = force_matrix[force_max_idx]

    
    power_list.append(2*np.sum(power_times_occ_matrix)/10e6)
    max_power_list.append(np.max(power_matrix)/10e6)
    force_list.append(2*np.sum(force_times_occ_matrix)/10e6)
    max_force_list.append(np.max(force_matrix)/10e6)
    
    
normalized_power_list = [100*(p-power_list[-1])/power_list[-1] for p in power_list]
normalized_max_power_list = [100*(p-max_power_list[-1])/max_power_list[-1] for p in max_power_list]
normalized_force_list = [100*(f-force_list[-1])/force_list[-1] for f in force_list]
normalized_max_force_list = [100*(f-max_force_list[-1])/max_force_list[-1] for f in max_force_list]


plt.figure(0)
plt.scatter(sorted_F5,normalized_power_list,marker='^',s=200,color='gold',label='mean power')
plt.scatter(sorted_F5,normalized_max_power_list,marker='^',s=75,color='orange',label='max power')
plt.scatter(sorted_F5,normalized_force_list,marker='*',s=200,color='orchid',label='mean force')
plt.scatter(sorted_F5,normalized_max_force_list,marker='*',s=75,color='purple',label='max force')
plt.xlabel('$|\widetilde{F_5}|$',fontsize=10)
plt.ylabel('Percent difference from shape C')
plt.legend()


plt.show()


