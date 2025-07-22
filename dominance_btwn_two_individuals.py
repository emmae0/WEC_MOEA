def dominance_btwn_two_individuals(param_1_1,param_1_2,param_2_1,param_2_2):
  if param_1_1 <= param_1_2 and param_2_1 <= param_2_2 and (param_1_1 < param_1_2 or param_2_1 < param_2_2):
    return 1
  elif param_1_1 >= param_1_2 and param_2_1 >= param_2_2 and (param_1_1 > param_1_2 or param_2_1 > param_2_2):
    return 2
  else:
    return 0
