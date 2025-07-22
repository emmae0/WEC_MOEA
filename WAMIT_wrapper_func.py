### WAMIT wrapper function: make sure it has converged
import csv
import numpy as np

def WAMIT_wrapper_func(geom_vector,make_mesh_figure,frequencies,NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55,beta_PTO,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup):

  ### look in lookup table to see if hydrodynamic coefficients saved for that shape. If it is, just output those. If it isn't, we will run the wrapper function and then save the coefficients in the lookuptable
  
  omega = frequencies[0]
  theta = betas[0]
  
  print ('starting')
  print ('geom vector:',geom_vector)
  print ('frequency:',omega)
  print ('theta:',theta)
    
  if lookup == 'yes':  
    
    if modes[0]=='1':
    
      lookup_filename = '/home/eedwards/MPS_WEC_opt/python_codes/lookup_table_w_X1_X3.csv'
      
    else:
      lookup_filename = '/home/eedwards/MPS_WEC_opt/python_codes/lookup_table.csv'

    theta_file = []  
    omega_file = []  
    geom_vector_file = []
    A55_file = []
    B55_file = []
    X5_file = []
    A11_file = []
    B11_file = []
    X1_file = []
    A33_file = []
    B33_file = []
    X3_file = []
    with open(lookup_filename) as csvfile:
      reader = csv.reader(csvfile,delimiter=' ')
      for row in reader:
        theta_file.append(row[0])
        omega_file.append(row[1])
        geom_vector_file.append(row[2])
        A55_file.append(row[3])
        B55_file.append(row[4])
        X5_file.append(row[5])  
        if modes[0]=='1':
          A11_file.append(row[6])
          B11_file.append(row[7])
          X1_file.append(row[8])
          A33_file.append(row[9])
          B33_file.append(row[10])
          X3_file.append(row[11])


    for file_idx in range(len(A55_file)):
      gv_file_split = str(geom_vector_file[file_idx]).split(',')
      this_geom_vector_file = []
      for gv_idx in range(len(geom_vector)):
        if gv_idx == 0:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][1:]))
        elif gv_idx == len(geom_vector)-1:
          this_geom_vector_file.append(float(gv_file_split[gv_idx][:-1]))
        else:
          this_geom_vector_file.append(float(gv_file_split[gv_idx]))
      count_same = 0
      for gv_idx in range(len(geom_vector)):
        if round(geom_vector[gv_idx],2) == round(this_geom_vector_file[gv_idx],2):
          count_same = count_same + 1
      if count_same == len(geom_vector) and omega == float(omega_file[file_idx]) and theta == float(theta_file[file_idx]):
        
        return (float(A55_file[file_idx]),float(B55_file[file_idx]),float(X5_file[file_idx]), float(A11_file[file_idx]),float(B11_file[file_idx]),float(X1_file[file_idx]),float(A33_file[file_idx]),float(B33_file[file_idx]),float(X3_file[file_idx]))
        
  ### get to this point if geom vector not in file. Run WAMIT wrapper
              
  # run with NN = NN_start
  NN = NN_start
  
  from check_num_panels import check_num_panels
    
  num_panels = check_num_panels(geom_vector,NN)
  
  print ('NN:',NN)
  print ('num panels:',num_panels)
    
  if num_panels > 38380:
    if lookup == 'yes':
      with open(lookup_filename,'a') as csvfile:
          writer = csv.writer(csvfile,delimiter=' ')
          if modes[0]=='1':
            writer.writerow([theta,omega,geom_vector,np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
          else:
            writer.writerow([theta,omega,geom_vector,np.nan, np.nan, np.nan])
  
    return np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
  
  from general_WAMIT_control import general_WAMIT_control
  xi_1_list,xi_3_list,xi_5_list,A11_list,A33_list,A55_list,B11_list,B33_list,B55_list,freq_list_out,X_1_list,X_3_list,X_5_list,mass = general_WAMIT_control(geom_vector,make_mesh_figure,frequencies,NN,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55,beta_PTO,k_PTO,ulen,IALTFRC)
  
  to_check_list = []
    
  # pitch  
  if modes[8]=='1':
    A55_1 = A55_list[-1]
    B55_1 = B55_list[-1]
    X5_1 = X_5_list[-1]
    to_check_list.append(A55_1)
    to_check_list.append(B55_1)
    to_check_list.append(X5_1)
  
  # run until converged
  while True:
  
    NN = NN + 1
    
    counter = 0
    
    from check_num_panels import check_num_panels
    
    num_panels = check_num_panels(geom_vector,NN)
    
    print ('NN:',NN)
    print ('num panels:',num_panels)
    
    if num_panels > 38380:
    
      if lookup == 'yes':
        with open(lookup_filename,'a') as csvfile:
          writer = csv.writer(csvfile,delimiter=' ')
          if modes[0]=='1':
            writer.writerow([theta,omega,geom_vector,np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
          else:
            writer.writerow([theta,omega,geom_vector,np.nan, np.nan, np.nan])
    
      return np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
  
    xi_1_list,xi_3_list,xi_5_list,A11_list,A33_list,A55_list,B11_list,B33_list,B55_list,freq_list_out,X_1_list,X_3_list,X_5_list,mass = general_WAMIT_control(geom_vector,make_mesh_figure,frequencies,NN,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55,beta_PTO,k_PTO,ulen,IALTFRC)

    to_check_list_2 = []
    # heave
    if modes[4]=='1':
      A33_2 = A33_list[0]
      B33_2 = B33_list[0]
      to_check_list_2.append(A33_2)
      to_check_list_2.append(B33_2)

    # surge
    if modes[0]=='1':
      A11_2 = A11_list[0]
      B11_2 = B11_list[0]
      to_check_list_2.append(A11_2)
      to_check_list_2.append(B11_2)
    
    # pitch  
    if modes[8]=='1':
      A55_2 = A55_list[-1]
      B55_2 = B55_list[-1]
      X5_2 = X_5_list[-1]
      to_check_list_2.append(A55_2)
      to_check_list_2.append(B55_2)
      to_check_list_2.append(X5_2)
      
    abs_error_list = [abs(to_check_list_2[idx]-to_check_list[idx]) for idx in range(len(to_check_list))]
    perc_error_list = [abs((to_check_list_2[idx]-to_check_list[idx])/to_check_list_2[idx]) for idx in range(len(to_check_list))]

    either_abs_or_perc_list = 0
    for idx in range(len(to_check_list)):
      #if abs_error_list[idx]< abs_error_wrapper or perc_error_list[idx]< perc_error_wrapper:
      if perc_error_list[idx]< perc_error_wrapper:
        either_abs_or_perc_list = either_abs_or_perc_list + 1
    
    if either_abs_or_perc_list == len(to_check_list):
      
      print ('breaking')
      print ('panel factor:',NN)
      print ('perc error:',perc_error_list)
      print ('geom vector:',geom_vector)
      #input()
      
      ### need to add to lookup table and then return values
      if lookup == 'yes':
        with open(lookup_filename,'a') as csvfile:
          writer = csv.writer(csvfile,delimiter=' ')
          if modes[0]=='1':
            writer.writerow([theta,omega,geom_vector,A55_list[0], B55_list[0], X_5_list[0], A11_list[0], B11_list[0], X_1_list[0], A33_list[0], B33_list[0], X_3_list[0]])
          else:
            writer.writerow([theta,omega,geom_vector,A55_list[0], B55_list[0], X_5_list[0]])


      return A55_list[0], B55_list[0], X_5_list[0], A11_list[0], B11_list[0], X_1_list[0], A33_list[0], B33_list[0], X_3_list[0]
      
    elif counter > 100:
    
      print ('counter > 100')
      sys.exit()
      
    else:
    
      print ('passing')
      print ('geom vector:',geom_vector)
      print ('panel factor:',NN)
      print ('perc error:',perc_error_list)
      print ('abs error:',abs_error_list)
      
      to_check_list = [to_check_list_2[idx] for idx in range(len(to_check_list_2))]

      NN = NN + 1
      counter = counter + 1
      pass
      
  
  
  
  
  
  
  
