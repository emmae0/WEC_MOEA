from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.integrate
import math

g = 9.81
pi = 3.14

def given_geom_vector_find_kW_F_and_I55(geom_vector,make_mesh_figure,omega,NN_start,modes,betas,water_depth,dir_name,rho,zCG,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,k):

  [s1,r2,r1,h,l,a2,b2,s2] = geom_vector

  ### find the center of buoyancy for the shape
  from calculate_xCB_and_zCB import calculate_xCB_and_zCB
  xCB,zCB = calculate_xCB_and_zCB(geom_vector)
  
  # set x center of gravity to be equal to x center of buoyancy
  xCG = xCB
  
  from calculate_surface_area import calculate_surface_area
  SW = calculate_surface_area(geom_vector)
  
  from calculate_volume import calculate_volume
  Volume = calculate_volume(geom_vector)
  
  # this code runs WAMIT
  from WAMIT_wrapper_func import WAMIT_wrapper_func
  I55_in = 100000 # note this isn't actually used
  beta_PTO_in = 100000 # this isn't actually used
  lookup='yes' # this says whether to use the lookup table (saved values of WAMIT outputs) or not. See 'MOEA_more_of_each' for the input
  A55_nd, B55_nd, X5_nd = WAMIT_wrapper_func(geom_vector,make_mesh_figure,[omega],NN_start,modes,betas,water_depth,dir_name,rho,xCG,zCG,I55_in,beta_PTO_in,k_PTO,perc_error_wrapper,abs_error_wrapper,ulen,IALTFRC,lookup)
  
  if str(A55_nd) == 'nan':
    return np.nan, np.nan, np.nan
  
  ### now calculate the quantities of intereste (F5, kW)
  
  # this part calculates the hydrostatic force
  S11 = l/3*((s1+r2)**3-(s1)**3)
  C55_1_nd = S11
  C55_2_nd = Volume*((zCB-s2)-(zCG-s2))
  C55 = rho*g*(C55_1_nd+C55_2_nd)
  
  # this dimensionalises the hydrodynamic coefficients
  A55 = A55_nd*rho*ulen**5
  B55 = B55_nd*rho*ulen**5*omega
  X5 = X5_nd*rho*g*ulen**3
  
  # calculte I55 to satisfy the resonance equation
  I55 = (C55+k_PTO)/omega**2-A55
  
  # set PTO damping coefficient equal to the radiation damping (to get maximum power)
  beta_PTO = B55
  
  # calculate motion
  xi_5_sq = X5**2/((C55+k_PTO-omega**2*(I55+A55))**2+omega**2*(beta_PTO+B55)**2)
  xi_5 = xi_5_sq**(1/2)
  
  # calculate force
  F5 = beta_PTO*omega*xi_5+k_PTO*xi_5

  # calculate power
  P_5 = 0.5*beta_PTO*omega**2*xi_5**2
  
  # calculate group velocity
  Vg = omega/k*0.5*(1+2*k*water_depth/(math.sinh(2*k*water_depth)))

  # incident power
  P_I = 1/2*rho*g*Vg
  
  # capture width
  W_5 = P_5/P_I
  
  # nondimensionalise
  kW = k*W_5
  klS = k*SW**(1/2)
  F5_nd = F5/(rho*omega**2/k**4)
  
  return kW,F5_nd,I55
