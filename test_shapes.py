#### calculate surface area from geom_vector
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

def test_shapes(geom_vector):

  s1,r2,r1,h,l,a2,b2,s2 = geom_vector
  
  ### define c1 and c2 curves
  x_c1_func = lambda s: (r1+s1)/2-a2+((r1-s1)/2)*(2*s-1)+a2*(2*(2*s-1)**2-1)
  x_c2_func = lambda s: (r1+s1+r2)/2-b2+((r1-s1-r2)/2)*(2*s-1)+b2*(2*(2*s-1)**2-1)
  
  s_list = np.linspace(0,1,1000)
  any_crosses = 0
  for s in s_list:
  
    c1_val = x_c1_func(s)
    c2_val = x_c2_func(s)
    c1_to_c2 = c2_val-c1_val
    if round(c1_to_c2,2) < 0:
      any_crosses = any_crosses + 1
  
  if any_crosses == 0:
    return 'yes'
    
  else:
    return 'no'



  
