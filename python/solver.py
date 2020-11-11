import numpy as np
# from rungekutta.rk4 import rk4
from finite_differences.finite_differences_o2 import *
import matplotlib.pyplot as plt

### constants
c = 1

### wave equation 0.2
# d phi / d t = phi
# d Pi / d t = c^2 d^2 phi / d x^2

### discretize space with uniform grid
Nt = 100
Nx = 100
deltax = 0.5
deltat = 0.1

def F(u,t): # time derivative
  # u is at one point in time here, hence only two-dimensional
  dudt = np.zeros([2,Nx])
  dudt[0,:] = u[1,:]

  # compute second derivative of phi along discretized x axis
  d2phidx2= np.zeros(Nx)
  for i in range(1,Nx-1): # weil nur die inneren Werte genommen werden sollen, die Randwerte = 0
    d2phidx2[i] = 1/deltax**2 * (u[0,i+1] - 2*u[0,i] + u[0,i-1])

  dudt[1,:] = c**2 * d2phidx2
  return dudt

def rk4(u0,deltat,T,F):
  phi0 = u0[0]
  Pi0 = u0[1]
  u = np.zeros([2, Nt, Nx])
  u[0,0,0] = phi0
  u[1,0,0] = Pi0
  t=1
  # time iteration
  for i in range(0,Nt-1):
    k1 = F(u[:,i,:], t)
    k2 = F(u[:,i,:] + 0.5*deltat* k1,t + 0.5*deltat)
    k3 = F(u[:,i,:] + 0.5*deltat* k2,t + 0.5*deltat)
    k4 = F(u[:,i,:] + deltat* k3,t + deltat)

    u[:,i+1,:] = u[:,i,:] + deltat*(1/6*k1 + 1/3*k2 +1/3*k3 + 1/6*k4)
  return u

u0 = [1,2]
u_rk4 = rk4(u0,deltat,Nt*deltat,F)
phi_t = u_rk4[0,:,5]

### Plotting the solutions in a position-time diagram ###
times = np.arange(0,Nt*deltat,deltat)
plt.plot(times, phi_t)
plt.xlabel('t')
plt.ylabel('phi(t)')
# # plt.axis('square')
# # plt.xlim(-12,12)
# # plt.ylim(49.94,50.008)
# plt.legend()
plt.grid(color = 'gainsboro')
plt.savefig("plot-phi_t.pdf")