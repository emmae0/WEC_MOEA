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

def make_gdf(geom_vector,NN,make_mesh_figure,rho,ulen):
  s1,r2,r1,h,l,a2,b2,s2 = geom_vector
  
  max_dim = max(s1+r2,h+abs(s2),l/2)

  # Header
  Header = 'WAMIT GDF file for MPS-UoP WEC optimization project'
  # Characteristic length scale, gravitational constant
  L1 = '%s %s'%(str(ulen),str(g))
  # Are x=0 and y=0 planes of symmetry?
  L2 = '0 0'

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

  Coord_str = ''


  ### formulate x and z lists for nodes: from top right, across the row and then to next row
  z_list_side = []
  x_list_side = []
  for NN_c_idx in range(NN_c):
    for NN_top_side_idx in range(NN_top_side):
      z_list_side.append(-NN_c_idx*h/NP_c-s2)
      x_list_side.append(x_c2_func(NN_c_idx/NP_c)-NN_top_side_idx*(x_c2_func(NN_c_idx/NP_c)-x_c1_func(NN_c_idx/NP_c))/NP_top_side)

  N_panels_1 = NP_c*NP_top_side
   
  if make_mesh_figure == 'yes':
    fig = plt.figure(2)
    ax = fig.add_subplot(projection='3d')

  ### Now put those into into the right order for WAMIT gdf file. start with panel in top right, work left along row and then start on next row
  for NP_c_idx in range(NP_c):
    for NP_top_side_idx in range(NP_top_side):
      x1 = x_list_side[NP_c_idx*NN_top_side+NP_top_side_idx]
      x2 = x_list_side[NP_c_idx*NN_top_side+NP_top_side_idx+1]
      x3 = x_list_side[NP_c_idx*NN_top_side+NN_top_side+NP_top_side_idx+1]
      x4 = x_list_side[NP_c_idx*NN_top_side+NN_top_side+NP_top_side_idx]
    
      z1 = z_list_side[NP_c_idx*NN_top_side+NP_top_side_idx]
      z2 = z_list_side[NP_c_idx*NN_top_side+NP_top_side_idx+1]
      z3 = z_list_side[NP_c_idx*NN_top_side+NN_top_side+NP_top_side_idx+1]
      z4 = z_list_side[NP_c_idx*NN_top_side+NN_top_side+NP_top_side_idx]
    
      y1 = -l/2
      y2 = -l/2
      y3 = -l/2
      y4 = -l/2
    
      ### makes a figure to visually represent the mesh
      if make_mesh_figure == 'yes':
        plt.figure(0)
        plt.plot([x1,x2,x3,x4,x1],[z1,z2,z3,z4,z1],color='b')
        plt.scatter([x1,x2,x3,x4],[z1,z2,z3,z4],color='b')
    
        ax.scatter(x1,y1,z1,color='k')
        ax.scatter(x2,y2,z2,color='k')
        ax.scatter(x3,y3,z3,color='k')
        ax.scatter(x4,y4,z4,color='k')
    
        ax.plot([x2,x1],[y2,y1],[z2,z1],color='k')
        ax.plot([x3,x2],[y3,y2],[z3,z2],color='k')
        ax.plot([x4,x3],[y4,y3],[z4,z3],color='k')
        ax.plot([x1,x4],[y1,y4],[z1,z4],color='k')
    
      new_str = ' '.join([str(x1),str(y1),str(z1),str(x2 ),str(y2),str(z2),str(x3),str(y3),str(z3),str(x4),str(y4),str(z4),'\n'])
    
      Coord_str = Coord_str + new_str
    
  ########### FACE 2 ##################

  l_top_front = l
  NP_top_front = math.ceil(l_top_front/l_panel)
  NN_top_front = NP_top_front+1

  N_panels_2 = NP_c*NP_top_front

  x_list_front = []
  y_list_front = []
  z_list_front = []

  for NN_c_idx in range(NN_c):
    for NN_top_front_idx in range(NN_top_front):
      y_list_front.append((-l/2)+(NN_top_front_idx)*(l)/NP_top_front)
      x_list_front.append(x_c1_func(NN_c_idx/NP_c))
      z_list_front.append(-NN_c_idx*h/NP_c-s2)
     
  
  ### Now put those into into the right order for WAMIT gdf file. start with panel in top right, work left along row and then start on next row
  for NP_c_idx in range(NP_c):
    for NP_top_front_idx in range(NP_top_front):
      x1 = x_list_front[NP_c_idx*NN_top_front+NP_top_front_idx]
      x2 = x_list_front[NP_c_idx*NN_top_front+NP_top_front_idx+1]
      x3 = x_list_front[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx+1]
      x4 = x_list_front[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx]
 
      z1 = z_list_front[NP_c_idx*NN_top_front+NP_top_front_idx]
      z2 = z_list_front[NP_c_idx*NN_top_front+NP_top_front_idx+1]
      z3 = z_list_front[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx+1]
      z4 = z_list_front[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx]
      
      
    
      y1 = y_list_front[NP_c_idx*NN_top_front+NP_top_front_idx]
      y2 = y_list_front[NP_c_idx*NN_top_front+NP_top_front_idx+1]
      y3 = y_list_front[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx+1]
      y4 = y_list_front[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx]
      

      if make_mesh_figure == 'yes':
        ax.scatter(x1,y1,z1,color='k')
        ax.scatter(x2,y2,z2,color='k')
        ax.scatter(x3,y3,z3,color='k')
        ax.scatter(x4,y4,z4,color='k')
        
       
    
        ax.plot([x2,x1],[y2,y1],[z2,z1],color='k')
        ax.plot([x3,x2],[y3,y2],[z3,z2],color='k')
        ax.plot([x4,x3],[y4,y3],[z4,z3],color='k')
        ax.plot([x1,x4],[y1,y4],[z1,z4],color='k')
        
    
      new_str = ' '.join([str(x1),str(y1),str(z1),str(x2 ),str(y2),str(z2),str(x3),str(y3),str(z3),str(x4),str(y4),str(z4),'\n'])
    
      Coord_str = Coord_str + new_str
      

  ########### FACE 3 #################
  N_panels_3 = NP_c*NP_top_side

  x_list_pos_side = []
  for NN_c_idx in range(NN_c):
    for NN_top_side_idx in range(NN_top_side):
      x_list_pos_side.append(x_c1_func(NN_c_idx/NP_c)+NN_top_side_idx*(x_c2_func(NN_c_idx/NP_c)-x_c1_func(NN_c_idx/NP_c))/NP_top_side)

  for NP_c_idx in range(NP_c):
    for NP_top_side_idx in range(NP_top_side):
      x1 = x_list_pos_side[NP_c_idx*NN_top_side+NP_top_side_idx]
      x2 = x_list_pos_side[NP_c_idx*NN_top_side+NP_top_side_idx+1]
      x3 = x_list_pos_side[NP_c_idx*NN_top_side+NN_top_side+NP_top_side_idx+1]
      x4 = x_list_pos_side[NP_c_idx*NN_top_side+NN_top_side+NP_top_side_idx]
    
      z1 = z_list_side[NP_c_idx*NN_top_side+NP_top_side_idx]
      z2 = z_list_side[NP_c_idx*NN_top_side+NP_top_side_idx+1]
      z3 = z_list_side[NP_c_idx*NN_top_side+NN_top_side+NP_top_side_idx+1]
      z4 = z_list_side[NP_c_idx*NN_top_side+NN_top_side+NP_top_side_idx]
    
      y1 = l/2
      y2 = l/2
      y3 = l/2
      y4 = l/2
    
      if make_mesh_figure == 'yes':
        ax.scatter(x1,y1,z1,color='k')
        ax.scatter(x2,y2,z2,color='k')
        ax.scatter(x3,y3,z3,color='k')
        ax.scatter(x4,y4,z4,color='k')
    
        ax.plot([x2,x1],[y2,y1],[z2,z1],color='k')
        ax.plot([x3,x2],[y3,y2],[z3,z2],color='k')
        ax.plot([x4,x3],[y4,y3],[z4,z3],color='k')
        ax.plot([x1,x4],[y1,y4],[z1,z4],color='k')
    
      new_str = ' '.join([str(x1),str(y1),str(z1),str(x2 ),str(y2),str(z2),str(x3),str(y3),str(z3),str(x4),str(y4),str(z4),'\n'])
    
      Coord_str = Coord_str + new_str
    
  ########### FACE 4 #####################   
  x_list_back = []
  y_list_back = []
  z_list_back = []

  N_panels_4 = NP_c*NP_top_front

  for NN_c_idx in range(NN_c):
    for NN_top_front_idx in range(NN_top_front):
      y_list_back.append(-l/2+(NP_top_front-NN_top_front_idx)*(l)/NP_top_front)
      x_list_back.append(x_c2_func(NN_c_idx/NP_c))
      z_list_back.append(-NN_c_idx*h/NP_c-s2)
 

  for NP_c_idx in range(NP_c):
    for NP_top_front_idx in range(NP_top_front):
      x1 = x_list_back[NP_c_idx*NN_top_front+NP_top_front_idx]
      x2 = x_list_back[NP_c_idx*NN_top_front+NP_top_front_idx+1]
      x3 = x_list_back[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx+1]
      x4 = x_list_back[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx]
 
      z1 = z_list_back[NP_c_idx*NN_top_front+NP_top_front_idx]
      z2 = z_list_back[NP_c_idx*NN_top_front+NP_top_front_idx+1]
      z3 = z_list_back[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx+1]
      z4 = z_list_back[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx]
    
      y1 = y_list_back[NP_c_idx*NN_top_front+NP_top_front_idx]
      y2 = y_list_back[NP_c_idx*NN_top_front+NP_top_front_idx+1]
      y3 = y_list_back[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx+1]
      y4 = y_list_back[NP_c_idx*NN_top_front+NN_top_front+NP_top_front_idx]
    
      if make_mesh_figure == 'yes':
        ax.scatter(x1,y1,z1,color='k')
        ax.scatter(x2,y2,z2,color='k')
        ax.scatter(x3,y3,z3,color='k')
        ax.scatter(x4,y4,z4,color='k')
    
        ax.plot([x2,x1],[y2,y1],[z2,z1],color='k')
        ax.plot([x3,x2],[y3,y2],[z3,z2],color='k')
        ax.plot([x4,x3],[y4,y3],[z4,z3],color='k')
        ax.plot([x1,x4],[y1,y4],[z1,z4],color='k')
    
        plt.figure(1)
        plt.plot([y1,y2,y3,y4,y1],[z1,z2,z3,z4,z1],color='b')
        plt.scatter([y1,y2,y3,y4],[z1,z2,z3,z4],color='b')
    
      new_str = ' '.join([str(x1),str(y1),str(z1),str(x2 ),str(y2),str(z2),str(x3),str(y3),str(z3),str(x4),str(y4),str(z4),'\n'])
    
      Coord_str = Coord_str + new_str

  
  ########### FACE 5 ##################
  x_list_roof = []
  y_list_roof = []

  N_panels_5 = NP_top_side*NP_top_front

  for NN_top_side_idx in range(NN_top_side):
    for NN_top_front_idx in range(NN_top_front):
      x_list_roof.append(s1+(NN_top_side_idx/NP_top_side)*r2)
      y_list_roof.append((NP_top_front-NN_top_front_idx)*(l)/NP_top_front-l/2)
 
  for NP_top_side_idx in range(NP_top_side):
    for NP_top_front_idx in range(NP_top_front):
      x1 = x_list_roof[NP_top_side_idx*NN_top_front+NP_top_front_idx]
      x2 = x_list_roof[NP_top_side_idx*NN_top_front+NP_top_front_idx+1]
      x3 = x_list_roof[NP_top_side_idx*NN_top_front+NN_top_front+NP_top_front_idx+1]
      x4 = x_list_roof[NP_top_side_idx*NN_top_front+NN_top_front+NP_top_front_idx]
 
      y1 = y_list_roof[NP_top_side_idx*NN_top_front+NP_top_front_idx]
      y2 = y_list_roof[NP_top_side_idx*NN_top_front+NP_top_front_idx+1]
      y3 = y_list_roof[NP_top_side_idx*NN_top_front+NN_top_front+NP_top_front_idx+1]
      y4 = y_list_roof[NP_top_side_idx*NN_top_front+NN_top_front+NP_top_front_idx]
    
      z1 = -s2
      z2 = -s2
      z3 = -s2
      z4 = -s2
    
      if make_mesh_figure == 'yes':
        ax.scatter(x1,y1,z1,color='k')
        ax.scatter(x2,y2,z2,color='k')
        ax.scatter(x3,y3,z3,color='k')
        ax.scatter(x4,y4,z4,color='k')
    
        ax.plot([x2,x1],[y2,y1],[z2,z1],color='k')
        ax.plot([x3,x2],[y3,y2],[z3,z2],color='k')
        ax.plot([x4,x3],[y4,y3],[z4,z3],color='k')
        ax.plot([x1,x4],[y1,y4],[z1,z4],color='k')
        
      new_str = ' '.join([str(x1),str(y1),str(z1),str(x2 ),str(y2),str(z2),str(x3),str(y3),str(z3),str(x4),str(y4),str(z4),'\n'])
    
      Coord_str = Coord_str + new_str


  total_N_Panels = N_panels_1 + N_panels_2 + N_panels_3 + N_panels_4 + N_panels_5 
  

  gdf_string = "\n".join([Header,L1,L2,str(total_N_Panels),Coord_str])

  if make_mesh_figure == 'yes':
    plt.figure(0)
    plt.xlabel('x')
    plt.ylabel('z')
    plt.xlim([-0.5,1.1*max_dim])
    plt.ylim([-1.1*max_dim,0])

    plt.figure(1)
    plt.xlabel('y')
    plt.ylabel('z')
    plt.xlim([-1.1*max_dim,1.1*max_dim])
    plt.ylim([-1.1*max_dim,0])

    plt.figure(2)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(0,1.1*max_dim)
    ax.set_ylim(-1.1*max_dim,1.1*max_dim)
    ax.set_zlim(-1.1*max_dim,0)

    ax.view_init(10,-45)
    
    plt.show()

  return gdf_string
  



