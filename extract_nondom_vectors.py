import numpy as np
import math

def extract_nondom_vectors(pop_V_vectors,pop_param_1,pop_param_2):

  from dominance_btwn_two_individuals import dominance_btwn_two_individuals
  
  #to keep track of dominated solutions
  dominated_solutions = []
  for dom_idx_1 in range(len(pop_V_vectors)):
    for dom_idx_2 in [x for x in range(len(pop_V_vectors)) if x != dom_idx_1]:
      # test dominance
      dom_result = dominance_btwn_two_individuals(pop_param_1[dom_idx_1],pop_param_1[dom_idx_2],pop_param_2[dom_idx_1],pop_param_2[dom_idx_2])
      # if 1 dominates 2, add 2 to dominated solutions
      if dom_result == 1:
        dominated_solutions.append(pop_V_vectors[dom_idx_2])
      # if 2 dominates 1, add 1 to dominated solutions
      elif dom_result == 2:
        dominated_solutions.append(pop_V_vectors[dom_idx_1])
  
  # Form Pareto Front by adding to current_PF if individual is NOT in dominated solutions
  PF_V_vectors = []
  PF_param_1 = []
  PF_param_2 = []
  for PF_idx in range(len(pop_V_vectors)):
    if pop_V_vectors[PF_idx] not in dominated_solutions:
      PF_V_vectors.append(pop_V_vectors[PF_idx])
      PF_param_1.append(pop_param_1[PF_idx])
      PF_param_2.append(pop_param_2[PF_idx])
      
  return PF_V_vectors,PF_param_1,PF_param_2
