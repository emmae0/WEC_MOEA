### make 2D figure of WEC
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

def plot_2D(geom_vector):

  # geometry vector that defines the shape -- see paper for definitions
  s1,r2,r1,h,l,a2,b2,s2 = geom_vector
  
  ### define c1 and c2 curves -- these are Chebyshev polynomials with coefficients a2 and b2. s goes from 0 to 1 along each curve.
  x_c1_func = lambda s: (r1+s1)/2-a2+((r1-s1)/2)*(2*s-1)+a2*(2*(2*s-1)**2-1)
  x_c2_func = lambda s: (r1+s1+r2)/2-b2+((r1-s1-r2)/2)*(2*s-1)+b2*(2*(2*s-1)**2-1)
  
  ### this turns it into x and z lists, which are what are output
  x_c1_list = []
  z_c1_list = []
  x_c2_list = []
  z_c2_list = []
  s_list = np.linspace(0,1,1000)
  for s in s_list:
    x_c1_list.append(x_c1_func(s))
    z_c1_list.append(-h*s)
    x_c2_list.append(x_c2_func(s))
    z_c2_list.append(-h*s)
    
  return x_c1_list,z_c1_list,x_c2_list,z_c2_list
    
  
