import numpy as np
import random

def select_parent_1(input_pop_V_vectors,input_pop_F5,input_pop_one_over_kW):
  from dominance_btwn_two_individuals import dominance_btwn_two_individuals
  
  # randomly select two individuals from the current population
  
  # select index of the vectors randomly
  tournament_vector_idx = np.random.randint(0,len(input_pop_V_vectors),size=2)
  # make list of vectors, 1/kW and F5 for chosen vectors
  tournament_vectors = [input_pop_V_vectors[idx] for idx in tournament_vector_idx]
  tournament_F5 = [input_pop_F5[idx] for idx in tournament_vector_idx]
  tournament_one_over_kW = [input_pop_one_over_kW[idx] for idx in tournament_vector_idx]
  
  # evaluate dominance 
  dom_result = dominance_btwn_two_individuals(tournament_F5[0],tournament_F5[1],tournament_one_over_kW[0],tournament_one_over_kW[1])
  
  if dom_result == 0:
    # if neither individual dominates, randomly choose one
    dom_result = random.randint(1,2)
  
  parent_V_vector = tournament_vectors[dom_result-1]
  parent_F5 = tournament_F5[dom_result-1]
  parent_one_over_kW = tournament_one_over_kW[dom_result-1]
  # -1 is there because dom result is 1 if frst one dominates (index 0)
    
  return parent_V_vector,parent_F5,parent_one_over_kW
