#### multi-objective optimization for MPS WEC
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math
import random
from given_geom_vector_find_kW_F_and_I55 import given_geom_vector_find_kW_F_and_I55


def MOEA(make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k,popsize,mut_prob,num_gens,possible_s1,possible_r2,possible_r1,possible_h,possible_a2,possible_b2,possible_s2,possible_l,show_opt_plots):

  #### if the list only has one entry, we aren't optimizing that parameter, so instead just set the value equal to the entry
  if len(possible_s1) == 1:
    s1 = possible_s1[0]
  if len(possible_r2) == 1:
    r2 = possible_r2[0]
  if len(possible_r1) == 1:
    r1 = possible_r1[0]
  if len(possible_h) == 1:
    h = possible_h[0]
  if len(possible_a2) == 1:
    a2 = possible_a2[0]
  if len(possible_b2) == 1:
    b2 = possible_b2[0]
  if len(possible_s2) == 1:
    s2 = possible_s2[0]
  if len(possible_l) == 1:
    l = possible_l[0]
    
  # tells us how many of each parameter we are considering
  num_each_param = [len(possible_s1),len(possible_r2),len(possible_r1),len(possible_h),len(possible_l),len(possible_a2),len(possible_b2),len(possible_s2)]
    
  ### define initial population
  
  # vector to keep track of what vectors have been tried
  V_vectors_tried = []
  
  init_pop_V_vectors = []
  for init_pop_idx in range(popsize):
    while len(init_pop_V_vectors) < popsize:
      # building partiulcar vector
      this_vector = []
      for init_pop_param_idx in range(len(num_each_param)):
        random_param_idx = random.randint(0,num_each_param[init_pop_param_idx]-1)
        this_vector.append(random_param_idx)
      
      # add if it's not in initial population yet
      if this_vector not in init_pop_V_vectors:
        init_pop_V_vectors.append(this_vector)
        
  # add to V_vectors_tried
  for V_idx in range(len(init_pop_V_vectors)):
    V_vectors_tried.append(init_pop_V_vectors[V_idx])
    
  ### get objective functions for each indivdual in initial population
  
  current_pop_V_vectors = []
  current_pop_one_over_kW = []
  current_pop_F5 = []
  current_pop_geom_vectors = []
  for vec_idx in range(len(init_pop_V_vectors)):
  
    print ('##################### initial population:',vec_idx,'#########################')
  
    # define geom vector
    
    geom_vector = [possible_s1[init_pop_V_vectors[vec_idx][0]],possible_r2[init_pop_V_vectors[vec_idx][1]],possible_r1[init_pop_V_vectors[vec_idx][2]],possible_h[init_pop_V_vectors[vec_idx][3]],possible_l[init_pop_V_vectors[vec_idx][4]],possible_a2[init_pop_V_vectors[vec_idx][5]],possible_b2[init_pop_V_vectors[vec_idx][6]],possible_s2[init_pop_V_vectors[vec_idx][7]]]

    ### test shape to make sure they don't cross over.
    from test_shapes import test_shapes
    shape_ok = test_shapes(geom_vector)
    
    if shape_ok == 'yes':
    
      try:

        kW,F5_nd,I55 = given_geom_vector_find_kW_F_and_I55(geom_vector,make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k)
        
        if str(kW) != 'nan' and I55 > 0:
          
          V_vector = init_pop_V_vectors[vec_idx]
    
          current_pop_V_vectors.append(V_vector)
          current_pop_one_over_kW.append(1/kW)
          current_pop_F5.append(F5_nd)
          current_pop_geom_vectors.append(geom_vector)
          

      except FileNotFoundError:
      
        geom_vector[1] = geom_vector[1]+0.01
        
        print (geom_vector)
        
        kW,F5_nd,I55 = given_geom_vector_find_kW_F_and_I55(geom_vector,make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k)
   
        V_vector = init_pop_V_vectors[vec_idx]
    
        current_pop_V_vectors.append(V_vector)
        current_pop_one_over_kW.append(1/kW)
        current_pop_F5.append(F5_nd)
        current_pop_geom_vectors.append(geom_vector)
    
  if show_opt_plots == 'yes':
    plt.figure(0)
    plt.scatter(current_pop_F5,[1/w for w in current_pop_one_over_kW],color='b')
    plt.xlabel('$F_5$')
    plt.ylabel('$kW$')

    plt.figure(1)
    plt.scatter(current_pop_F5,current_pop_one_over_kW,color='b')
    plt.xlabel('$F_5$')
    plt.ylabel('$1/kW$')

  ### extract nondominated individuals from the initial population to form the initial Pareto Front
  
  from extract_nondom_vectors import extract_nondom_vectors
  
  current_PF_V_vectors,current_PF_one_over_kW,current_PF_F5 = extract_nondom_vectors(current_pop_V_vectors,current_pop_one_over_kW,current_pop_F5)
  
  if show_opt_plots == 'yes':
    plt.figure(0)
    plt.scatter(current_PF_F5,[1/w for w in current_PF_one_over_kW],color='g')

    plt.figure(1)
    plt.scatter(current_PF_F5,current_PF_one_over_kW,color='g')  
    
  # get geom_vectors for PF vectors
  
  current_PF_geom_vectors = []
  
  for vec_idx in range(len(current_PF_V_vectors)):
    geom_vector = [possible_s1[current_PF_V_vectors[vec_idx][0]],possible_r2[current_PF_V_vectors[vec_idx][1]],possible_r1[current_PF_V_vectors[vec_idx][2]],possible_h[current_PF_V_vectors[vec_idx][3]],possible_l[current_PF_V_vectors[vec_idx][4]],possible_a2[current_PF_V_vectors[vec_idx][5]],possible_b2[current_PF_V_vectors[vec_idx][6]],possible_s2[current_PF_V_vectors[vec_idx][7]]]
    current_PF_geom_vectors.append(geom_vector)
    
  ### start algorithm loop
  
  for loop_idx in range(num_gens):
    
    print ('##################### generation:',loop_idx,'#########################')
  
    # if we've tried all possible vectors, return population etc now
    
    num_vectors_poss = 1
    for param_idx in range(len(num_each_param)):
      num_vectors_poss = num_vectors_poss*num_each_param[param_idx]
      
    if len(V_vectors_tried) == num_vectors_poss:
      return current_pop_V_vectors,current_pop_F5,current_pop_one_over_kW,current_pop_geom_vectors,current_PF_V_vectors,current_PF_F5,current_PF_one_over_kW,current_PF_geom_vectors
      
    ### one parent is chosen from current population
    
    from select_parent_1 import select_parent_1
    
    first_parent_V_vector,first_parent_F5,first_parent_one_over_kW = select_parent_1(current_pop_V_vectors,current_pop_F5,current_pop_one_over_kW)
    
    
    ### other parent is chosen from current PF
    
    second_parent_idx = random.randint(0,len(current_PF_V_vectors)-1)
    second_parent_V_vector,second_parent_F5,second_parent_one_over_kW = current_PF_V_vectors[second_parent_idx],current_PF_F5[second_parent_idx],current_PF_one_over_kW[second_parent_idx]
    
    ### child is created using parents
    
    from create_offspring import create_offspring
    
    child_V_vector,child_F5,child_one_over_kW,child_geom_vector,V_vectors_tried = create_offspring(first_parent_V_vector,first_parent_F5,first_parent_one_over_kW,second_parent_V_vector,second_parent_F5,second_parent_one_over_kW,make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k,mut_prob,possible_s1,possible_r2,possible_r1,possible_h,possible_a2,possible_b2,possible_s2,possible_l,V_vectors_tried)

    ### update population
    
    if child_V_vector != 'no':
    
      current_pop_V_vectors.append(child_V_vector)
      current_pop_F5.append(child_F5)
      current_pop_one_over_kW.append(child_one_over_kW)
      current_pop_geom_vectors.append(child_geom_vector)
    
    ### update PF
    
    current_PF_V_vectors,current_PF_one_over_kW,current_PF_F5 = extract_nondom_vectors(current_pop_V_vectors,current_pop_one_over_kW,current_pop_F5)

    current_PF_geom_vectors = []
  
    for vec_idx in range(len(current_PF_V_vectors)):
      geom_vector = [possible_s1[current_PF_V_vectors[vec_idx][0]],possible_r2[current_PF_V_vectors[vec_idx][1]],possible_r1[current_PF_V_vectors[vec_idx][2]],possible_h[current_PF_V_vectors[vec_idx][3]],possible_l[current_PF_V_vectors[vec_idx][4]],possible_a2[current_PF_V_vectors[vec_idx][5]],possible_b2[current_PF_V_vectors[vec_idx][6]],possible_s2[current_PF_V_vectors[vec_idx][7]]]
      current_PF_geom_vectors.append(geom_vector)
    
  return current_pop_V_vectors,current_pop_F5,current_pop_one_over_kW,current_pop_geom_vectors,current_PF_V_vectors,current_PF_F5,current_PF_one_over_kW,current_PF_geom_vectors
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
