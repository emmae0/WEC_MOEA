### python code to define the geometry of the WEC and make the mesh
from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

### fixed parameters
pi = 3.14
g = 9.81

def check_num_panels(geom_vector,NN):
  s1,r2,r1,h,l,a2,b2,s2 = geom_vector
  
  max_dim = max(s1+r2,h+abs(s2),l/2)


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
  # length of top is just r2
  l_top = r2

  ######## FACE 1 ############

  ### determine how many nodes and panel on each side
  # calculate side with min length
  min_l = min(l_top,l_c1,l_c2)
  # there are NN nodes on that smallest side. determine length of panel
  l_panel = min_l/NN
  # now determine how many nodes and panels on each side
  NP_top_side = math.ceil(l_top/l_panel)
  NN_top_side = NP_top_side+1
  NP_c = math.ceil(min(l_c1,l_c2)/l_panel)
  NN_c = NP_c+1



  N_panels_1 = NP_c*NP_top_side
   
  l_top_front = l
  NP_top_front = math.ceil(l_top_front/l_panel)
  NN_top_front = NP_top_front+1

  N_panels_2 = NP_c*NP_top_front



  N_panels_3 = NP_c*NP_top_side


  N_panels_4 = NP_c*NP_top_front


  #N_panels_5 = NP_top_side*NP_top_front



  total_N_Panels = N_panels_1 + N_panels_2 + N_panels_3 + N_panels_4 #+ N_panels_5 
  

  return total_N_Panels
  






