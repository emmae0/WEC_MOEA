#### calculate surface area from geom_vector
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

def calculate_surface_area(geom_vector):

  s1,r2,r1,h,l,a2,b2,s2 = geom_vector
  
  ### define c1 and c2 curves
  x_c1_func = lambda s: (r1+s1)/2-a2+((r1-s1)/2)*(2*s-1)+a2*(2*(2*s-1)**2-1)
  x_c2_func = lambda s: (r1+s1+r2)/2-b2+((r1-s1-r2)/2)*(2*s-1)+b2*(2*(2*s-1)**2-1)
  
  ### derivative functions (for calculating arclength, volume, etc etc)
  x_c1_prime_func = lambda s: r1-s1+4*a2*s
  x_c2_prime_func = lambda s: r1-s1-r2+4*b2*s
  
  ### calculate arclengths for c1 and c2
  # s list for derivative calculation
  s_list = np.linspace(0,1,100)
  # inside deriv lists
  l_c1_inside = [(h**2+(x_c1_prime_func(s))**2)**(1/2) for s in s_list]
  l_c2_inside = [(h**2+(x_c2_prime_func(s))**2)**(1/2) for s in s_list]
  # calculate arclength by integrating above lists
  l_c1 = scipy.integrate.trapz(l_c1_inside,s_list)
  l_c2 = scipy.integrate.trapz(l_c2_inside,s_list)
  
  s_list = np.linspace(0,1,100)
  horiz_lines_list = []
  for s in s_list:
  
    c1_val = x_c1_func(s)
    c2_val = x_c2_func(s)
    c1_to_c2 = c2_val-c1_val
    horiz_lines_list.append(c1_to_c2)
    
  area_2d = scipy.integrate.trapz(horiz_lines_list,s_list)*h
  
  surface_area = 2*area_2d+(l_c1+l_c2)*l
  
  return surface_area



  
