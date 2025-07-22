#### calculate surface area from geom_vector
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

def calculate_xCB_and_zCB(geom_vector):

  s1,r2,r1,h,l,a2,b2,s2 = geom_vector
  
  ### define c1 and c2 curves
  x_c1_func = lambda s: (r1+s1)/2-a2+((r1-s1)/2)*(2*s-1)+a2*(2*(2*s-1)**2-1)
  x_c2_func = lambda s: (r1+s1+r2)/2-b2+((r1-s1-r2)/2)*(2*s-1)+b2*(2*(2*s-1)**2-1)
  
  ### derivative functions (for calculating arclength, volume, etc etc)
  x_c1_prime_func = lambda s: r1-s1+4*a2*s
  x_c2_prime_func = lambda s: r1-s1-r2+4*b2*s
 
  
  s_list = np.linspace(0,1,1000)
  inside_sum_horiz_list = []
  horiz_lines_list = []
  inside_sum_vert_list = []
  for s in s_list:
  
    c1_val = x_c1_func(s)
    c2_val = x_c2_func(s)
    c1_to_c2 = c2_val-c1_val
    horiz_lines_list.append(c1_to_c2)
    inside_sum_horiz_list.append(c1_to_c2*((c1_val+c2_val)/2))
    inside_sum_vert_list.append(-s*h*c1_to_c2)
    

  area_2d = scipy.integrate.trapz(horiz_lines_list,s_list)*h
  
  x_CB = scipy.integrate.trapz(inside_sum_horiz_list,s_list)*h/area_2d
  
  z_CB = scipy.integrate.trapz(inside_sum_vert_list,s_list)*h/area_2d
  

  return x_CB,z_CB




  
