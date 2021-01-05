import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

from finite_differences.example_functions import *
from results_plotting import *
### constants
c = 1

### wave equation 0.2
# d phi / d t = phi
# d Pi / d t = c^2 d^2 phi / d x^2

### discretize time and space with uniform grid
def gridmaker(endT,nt,endX,nx):
  dt = (endT/nt)
  timevalues = np.linspace(0,endT,nt)
  dx = (endX/nx)
  xvalues = np.linspace(0,endX,nx)
  # print('dt = %.3f, dx = %.3f' %(dt,dx))
  return dx,timevalues,dx, xvalues


def wave_evolution1D(phi0,Pi0,timevalues,xvalues,bc):
  Nt = len(timevalues)
  Nx = len(xvalues)
  deltax = xvalues[Nx-1]/Nx
  deltat = timevalues[Nt-1]/Nt
  phi = np.zeros([Nt,Nx+2])
  Pi = np.zeros([Nt,Nx+2])
  ### first time step
  phi[0,1:Nx+1] = phi0
  Pi[0,1:Nx+1] = Pi0

  def rhs(phi,Pi,t):
    if bc == "periodic":
      phi[0] = phi[-3]
      phi[-1] = phi[2]
      Pi[0] = Pi[-3]
      Pi[-1] = Pi[2]
    elif bc == "Dirichlet": # like a string
      # phi[0] = 0
      # phi[-1] = 0
      # phi[1] = 0
      # phi[-2] = 0
      # Pi[0] = 0
      # Pi[1] = 0
      # Pi[-2] = 0
      # Pi[-1] = 0
      phi[0] = - phi[2]
      phi[-1] = - phi[-3]
      Pi[0] = - Pi[2]
      Pi[-1] = - Pi[-3]
    elif bc == "vonNeumann": # like a reflected water wave
      phi[0] = phi[1]
      phi[-1] = phi[-2]
      Pi[0] = Pi[1]
      Pi[-1] = Pi[-2]
    elif bc == "open_i":     # variant (i)     WTF????
      phi[0] = 2 * phi[1] - phi[2]
      phi[-1] = 2 * phi[-2] - phi[-3]
      Pi[0] = 2 * Pi[1] - Pi[2]
      Pi[-1] = 2* Pi[-2] - Pi[-3]
    elif bc == "open_ii":         # variant (ii)
      phi[0] = phi[2] - 2*Pi[1] * deltax/c
      phi[-1] = phi[-3] - 2*Pi[-2] * deltax/c
      Pi[0] = Pi[2] - 2*phi[1] * deltax/c
      Pi[-1] = Pi[-3] - 2*phi[-2] * deltax/c
      pass

    # compute second spatial derivative (d^2 phi / dx^2) with FD
    d2phidx2= np.zeros(Nx+2)
    for ix in range(1,Nx+1): # computing only inner points
      d2phidx2[ix] = 1/deltax**2 * (phi[ix+1] - 2*phi[ix] + phi[ix-1])

    dphidt = Pi
    dPidt = c**2 * d2phidx2
    return dphidt, dPidt

  t = 1 #dummy value
  # time iteration (RK4 method)
  for i in range(0,Nt-1):
    k1_phi, k1_Pi  = rhs(phi[i,:], Pi[i], t)
    k2_phi, k2_Pi = rhs(phi[i,:] + 0.5*deltat*k1_phi,Pi[i,:] + 0.5*deltat*k1_Pi,t + 0.5*deltat)
    k3_phi, k3_Pi = rhs(phi[i,:] + 0.5*deltat*k2_phi,Pi[i,:] + 0.5*deltat*k2_Pi ,t + 0.5*deltat)
    k4_phi, k4_Pi = rhs(phi[i,:] + deltat*k3_phi,Pi[i,:] + deltat*k3_Pi ,t + deltat)

    phi[i+1,:] = phi[i,:] + deltat*(1/6*k1_phi + 1/3*k2_phi +1/3*k3_phi + 1/6*k4_phi)
    Pi[i+1,:] = Pi[i,:] + deltat*(1/6*k1_Pi + 1/3*k2_Pi +1/3*k3_Pi + 1/6*k4_Pi)
  return phi[:,1:Nx+1], Pi[:,1:Nx+1] # return only inner points


#--------------------- take a look at the energy ------
def energy(q,p):        #calculate energy from position q(phi) and inertia p(pi)
    m=1
    E = 0.5* p**2 / m + 0.5* q**2 *m
    return E
def total_energy(phi,pi):
    (rows,columns) = np.shape(phi)
    Etotal = np.zeros(rows)
    E = energy(phi,pi)
    for i in range(0,rows):   # for all times sum up individual energies
        #divide by number of columns to make E independent of Nx
        Etotal[i] = sum(E[i,1:Nx+1])/columns # do not consider ghost points
    Etotal = Etotal/max(Etotal)
    return Etotal

# -------------------- now, do it ---------------
if __name__ == "__main__":
    endT = 1
    Nt = 400
    endX = 1
    Nx = 200
    sigma = 0.005
    mu = 0.7
    width= 0.2
    k = 1

    deltat, timevalues, deltax, xvalues = gridmaker(endT,Nt,endX,Nx)
    # courant = c * deltat / deltax
    # print("courant number = %.2f" % courant)

    ### choose f_4, f_5, g_a (for latter specify a = ...) or gaussian here (for latter specify sigma and mu)

    # Phi0 = f_4(xvalues)
    # Pi0  = - f_4_prime(xvalues)
    # Phi0 = gaussian(xvalues,sigma,mu)
    # Pi0  = -gaussian_drv(xvalues,sigma,mu)
    # Phi0 = squares(xvalues, k)
    # Pi0  = -squares_drv(xvalues,k)
    Phi0 = f_triangle(xvalues,width/2,mu)
    Pi0 = -f_triangle_drv(xvalues,width/2,mu)
    # Phi0 = g_a(xvalues,20)#
    # Pi0 = 3*np.zeros(len(phi0))#- g_a_prime(xvalues,20)

    # print(Phi0)
    # print(Pi0)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues, "periodic")

    Etotal = total_energy(Phi,Pi)
    # Nt_plot = 7 # how many snap shots are plotted
    plot_energy_evolution(Etotal,timevalues)
    plot_xt_evolution_heatmap(timevalues,xvalues,Phi)
    # plot_animation(xvalues, timevalues, Phi, Pi,'gif')

    # save as csv file
    # np.savetxt("results.csv", Phi, delimiter = ',', fmt = '%.6e')
