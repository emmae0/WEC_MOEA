import random

# biased coin
def flip(p):
  return 'H' if random.random() < p else 'T'

def create_offspring(first_parent_V_vector,first_parent_F5,first_parent_one_over_kW,second_parent_V_vector,second_parent_F5,second_parent_one_over_kW,make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k,mut_prob,possible_s1,possible_r2,possible_r1,possible_h,possible_a2,possible_b2,possible_s2,possible_l,V_vectors_tried):

  lookup = 'yes'

  num_each_param = [len(possible_s1),len(possible_r2),len(possible_r1),len(possible_h),len(possible_l),len(possible_a2),len(possible_b2),len(possible_s2)]
  num_params = len(num_each_param)

  # make list with indices of parameters that have more than one possibility 
  non_zero_params = []
  for param_idx in range(num_params):
    if num_each_param[param_idx] > 1:
      non_zero_params.append(param_idx)
      

  ### First, do crossover, then flip coin and do mutation if heads. Then, see if V vector is in V vectors tried. If it is, try again, If it isn't, that's child
  
  break_sig = 0 # this becomes 1 when we have a V vector that hasn't been used
  
  counter = 0
  while break_sig < 1:
  
    # to make sure doesn't get caught in loop
    counter = counter + 1
    if counter > 1000:
      print ('counter > 1000')
      input()
      
    ### crossover  
      
    child_V_vector_co = [0]*num_params
      
    # if there is only one parameter that can be chagned, child equals parent 1
    if len(non_zero_params) == 1:
      child_V_vector_co[non_zero_params[0]] = first_parent_V_vector[non_zero_params[0]]
      
    else:
      # randomly choose which parameter to do crossover at (to the left, parent 1; to the right, parent 2)
      crossover_idx = random.choice(non_zero_params[:-1])
      for coidx in range(num_params):
        if coidx < crossover_idx + 1:
          child_V_vector_co[coidx] = first_parent_V_vector[coidx]
        else:
          child_V_vector_co[coidx] = second_parent_V_vector[coidx]
    
    ### mutation
    mut = flip(mut_prob)
    
    if mut == 'H':
      
      mutated_child = []
      
      mutation_idx = random.choice(non_zero_params)
      param_idx_mut = random.randint(0,num_each_param[mutation_idx]-1)
      
      mutated_child = child_V_vector_co
      mutated_child[mutation_idx] = param_idx_mut
      
      if mutated_child not in V_vectors_tried:
        child_V_vector = mutated_child
        break_sig = 1
      else:
        break_sig = 0 # try again
        
    else:
    
      if child_V_vector_co not in V_vectors_tried:
        child_V_vector = child_V_vector_co
        break_sig = 1
      
      else:
      
        # if crossover is in vectors tried, do mutation anyways (even if didn't flip)
        
        mutation_idx = random.choice(non_zero_params)
        
        # check if all values of the chosen parameter are already in V vectors tried
        yes_count_mut = 0
        for param_idx_mut_test in range(num_each_param[mutation_idx]):
          test_V = []
          test_V = child_V_vector_co
          test_V[mutation_idx] = param_idx_mut_test
          if test_V in V_vectors_tried:
            yes_count_mut = yes_count_mut + 1
            
        if yes_count_mut == num_each_param[mutation_idx]:
        
          mutated_child = [0]*num_params
          for each_idx in range(num_params):
            if each_idx not in non_zero_params:
              mutated_child[each_idx] = child_V_vector_co[each_idx]
            else:
              param_idx_mut_ev = random.randint(0,num_each_param[each_idx]-1)
              mutated_child[each_idx] = param_idx_mut_ev
              
        else:
        
          param_idx_mut = random.randint(0,num_each_param[mutation_idx]-1)
          
          mutated_child = []
          mutated_child = child_V_vector_co
          mutated_child[mutation_idx] = param_idx_mut
          
        if mutated_child not in V_vectors_tried:
          child_V_vector = mutated_child
          break_sig = 1
        else:
          break_sig = 0
          
  
  V_vectors_tried.append(child_V_vector)
  
  geom_vector = [possible_s1[child_V_vector[0]],possible_r2[child_V_vector[1]],possible_r1[child_V_vector[2]],possible_h[child_V_vector[3]],possible_l[child_V_vector[4]],possible_a2[child_V_vector[5]],possible_b2[child_V_vector[6]],possible_s2[child_V_vector[7]]]
  
  ### test shape to make sure shape doesn't cross itself
  from test_shapes import test_shapes
  shape_ok = test_shapes(geom_vector)
  
  if shape_ok == 'yes':
    from given_geom_vector_find_kW_F_and_I55 import given_geom_vector_find_kW_F_and_I55
    
    try:
  
      kW,F5_nd,I55 = given_geom_vector_find_kW_F_and_I55(geom_vector,make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k)
    
      if str(kW) != 'nan' and I55 > 0:
        print ('Child: all good')
        return child_V_vector,F5_nd,1/kW,geom_vector,V_vectors_tried
        
      else:
        return 'no','no','no','no',V_vectors_tried
    
    except FileNotFoundError:
    
      geom_vector[1] = geom_vector[1]+0.01
        
      kW,F5_nd,I55 = given_geom_vector_find_kW_F_and_I55(geom_vector,make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k)

      if str(kW) != 'nan' and I55> 0:    
    
        print ('Child: Added 0.01')
        print (child_V_vector,F5_nd,1/kW,geom_vector,V_vectors_tried)
        return child_V_vector,F5_nd,1/kW,geom_vector,V_vectors_tried
        
      else:
        print ('Child: Added 0.01 but did not work')
        return 'no','no','no','no',V_vectors_tried
      
    
  else:
    print ('Child: shape not ok')
    return 'no','no','no','no',V_vectors_tried
      

  
  
  











