import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import math

from finite_differences.example_functions import *
from results_plotting import *
### constants
c = 1

### wave equation 0.2
# d phi / d t = phi
# d Pi / d t = c^2 d^2 phi / d x^2

### discretize time and space with uniform grid
def gridmaker(endT,nt,endX,nx,startX=0):
  dt = (endT/nt)
  timevalues = np.linspace(0,endT,nt)
  dx = ((endX-startX)/nx)
  xvalues = np.linspace(startX,endX,nx)
  print('dt = %.3f, dx = %.3f' %(dt,dx))
  return dt,timevalues,dx,xvalues

def wave_evolution1D(phi0,Pi0,timevalues,xvalues,bc,potential):
  Nt = len(timevalues)
  Nx = len(xvalues)
  deltax = (xvalues[Nx-1]-xvalues[0])/Nx
  deltat = timevalues[Nt-1]/Nt
  phi = np.zeros([Nt,Nx+2])     # add 1 ghost point at each end
  Pi = np.zeros([Nt,Nx+2])
  ### first time step
  phi[0,1:Nx+1] = phi0          # fill inner points with given Phi0 and Pi0
  Pi[0,1:Nx+1] = Pi0
  potential = np.insert(potential,0,potential[0])   # expand potential to fit to ghosts
  potential = np.append(potential,potential[-1])
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
    elif bc == "open_i":     # variant (i)
      phi[0] = 2 * phi[1] - phi[2]
      phi[-1] = 2 * phi[-2] - phi[-3]
      Pi[0] = 2 * Pi[1] - Pi[2]
      Pi[-1] = 2* Pi[-2] - Pi[-3]
    elif bc == "open_ii":         # variant (ii)
      phi[0] = phi[2] - 2*Pi[1] * deltax/c
      phi[-1] = phi[-3] - 2*Pi[-2] * deltax/c
      Pi[0] = Pi[2] - 2*phi[1] * deltax/c
      Pi[-1] = Pi[-3] - 2*phi[-2] * deltax/c
    elif bc == "open_iii":         # variant (iii)
      phi[0] = - deltax*(2*Pi[1] - Pi[2])/c + phi[1]
      phi[-1] = - deltax*(2*Pi[-2] - Pi[-3])/c + phi[-2]
      Pi[0] = (phi[0] - phi[1])/deltax
      Pi[-1] = (phi[-1] - phi[-2])/deltax

    # compute second spatial derivative (d^2 phi / dx^2) with FD
    d2phidx2= np.zeros(Nx+2)
    for ix in range(1,Nx+1): # computing only inner points
      d2phidx2[ix] = 1/deltax**2 * (phi[ix+1] - 2*phi[ix] + phi[ix-1])

    dphidt = Pi
    dPidt = c**2 * d2phidx2 + potential*phi
    # print(dPidt[300])
    # print(c**2 * d2phidx2[300])
    return dphidt, dPidt

  t = 1 #dummy value
  # time iteration (RK4 method)
  for i in range(0,Nt-1):
    k1_phi, k1_Pi = rhs(phi[i,:], Pi[i], t)
    k2_phi, k2_Pi = rhs(phi[i,:] + 0.5*deltat*k1_phi,Pi[i,:] + 0.5*deltat*k1_Pi,t + 0.5*deltat)
    k3_phi, k3_Pi = rhs(phi[i,:] + 0.5*deltat*k2_phi,Pi[i,:] + 0.5*deltat*k2_Pi ,t + 0.5*deltat)
    k4_phi, k4_Pi = rhs(phi[i,:] + deltat*k3_phi,Pi[i,:] + deltat*k3_Pi ,t + deltat)

    phi[i+1,:] = phi[i,:] + deltat*(1/6*k1_phi + 1/3*k2_phi +1/3*k3_phi + 1/6*k4_phi)
    Pi[i+1,:] = Pi[i,:] + deltat*(1/6*k1_Pi + 1/3*k2_Pi +1/3*k3_Pi + 1/6*k4_Pi)
  return phi[:,1:Nx+1], Pi[:,1:Nx+1] # return only inner points

def wave_evolution1D_4th_order(phi0,Pi0,timevalues,xvalues,bc,potential):
  Nt = len(timevalues)
  Nx = len(xvalues)
  deltax = xvalues[Nx-1]/Nx
  deltat = timevalues[Nt-1]/Nt
  phi = np.zeros([Nt,Nx+4])     # add 4 ghost point, 2 on each side of x domain
  Pi = np.zeros([Nt,Nx+4])       # add 4 ghost point, 2 on each side
  ### first time step
  phi[0,2:Nx+2] = phi0
  Pi[0,2:Nx+2] = Pi0
  potential = np.insert(potential,0,potential[0])
  potential = np.insert(potential,0,potential[1])
  potential = np.append(potential,potential[-1])
  potential = np.append(potential,potential[-2])
  def rhs(phi,Pi,t):
    if bc == "periodic":
      phi[0] = phi[-5]  # filling ghost points at the beginning
      phi[1] = phi[-4]
      phi[-1] = phi[4]  # filling ghost points at the end
      phi[-2] = phi[3]
      Pi[0] = Pi[-5]    # filling ghost points at the beginning
      Pi[1] = Pi[-4]
      Pi[-1] = Pi[4]    # filling ghost points at the end
      Pi[-2] = Pi[3]

    # compute second spatial derivative (d^2 phi / dx^2) with FiniteDiff
    d2phidx2= np.zeros(Nx+4)
    for ix in range(2,Nx+2): # computing only inner points
      d2phidx2[ix] = 1/(12*(deltax**2)) * (
      -phi[ix+2] + 16*phi[ix+1] -30*phi[ix]  +16*phi[ix-1] -phi[ix-2]
      )

    dphidt = Pi
    # print(np.shape(d2phidx2))
    # print(np.shape(phi),np.shape(potential))
    dPidt = c**2 * d2phidx2 - potential*phi
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
  return phi[:,2:Nx+2], Pi[:,2:Nx+2]    # return only inner points

def wave_evolution1D_6th_order(phi0,Pi0,timevalues,xvalues,bc,potential):
  Nt = len(timevalues)
  Nx = len(xvalues)
  deltax = xvalues[Nx-1]/Nx
  deltat = timevalues[Nt-1]/Nt
  phi = np.zeros([Nt,Nx+6])     # add 6 ghost point, 3 on each side of x domain
  Pi = np.zeros([Nt,Nx+6])       # add 6 ghost point, 3 on each side
  ### first time step
  phi[0,3:Nx+3] = phi0
  Pi[0,3:Nx+3] = Pi0
  potential = np.insert(potential,0,potential[0])
  potential = np.insert(potential,0,potential[1])
  potential = np.insert(potential,0,potential[2])
  potential = np.append(potential,potential[-1])
  potential = np.append(potential,potential[-2])
  potential = np.insert(potential,0,potential[-3])
  def rhs(phi,Pi,t):
    if bc == "periodic":
      phi[0] = phi[-7]  # filling ghost points at the beginning
      phi[1] = phi[-6]
      phi[2] = phi[-5]
      phi[-1] = phi[6]  # filling ghost points at the end
      phi[-2] = phi[5]
      phi[-3] = phi[4]
      Pi[0] = Pi[-7]    # filling ghost points at the beginning
      Pi[1] = Pi[-6]
      Pi[2] = Pi[-5]
      Pi[-1] = Pi[6]    # filling ghost points at the end
      Pi[-2] = Pi[5]
      Pi[-3] = Pi[4]

    # compute second spatial derivative (d^2 phi / dx^2) with FiniteDiff
    d2phidx2= np.zeros(Nx+6)
    for ix in range(3,Nx+3): # computing only inner points
      d2phidx2[ix] = 1/(90*(deltax**2)) * (
      phi[ix+3] -13.5*phi[ix+2] + 135*phi[ix+1] -245*phi[ix]
        +135*phi[ix-1] -13.5*phi[ix-2]+phi[ix+3]
      )

    dphidt = Pi
    # print(np.shape(d2phidx2))
    # print(np.shape(phi),np.shape(potential))
    dPidt = c**2 * d2phidx2 - potential*phi
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
  return phi[:,3:Nx+3], Pi[:,3:Nx+3] # return only inner points


#--------------------- take a look at the energy -------------------
#-------------------------------------------------------------------
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


# -------------------- little helper function ---------------
def IVmaker(func,xvalues,sigma,mu,width,k,ampl):
  funcDict = {"sine4":(f_4(xvalues),- f_4_prime(xvalues))
  ,"sine5":(f_5(xvalues),- f_5_prime(xvalues))
  ,"gauss":(gaussian(xvalues,sigma,mu,ampl),-gaussian_drv(xvalues,sigma,mu,ampl))
  ,"square":(squares(xvalues, k),-squares_drv(xvalues,k))
  # ,"triangle":(f_triangle(xvalues,width/2,mu),-f_triangle_drv(xvalues,width/2,mu))
  }
  return funcDict[func]


# ------------------- potential ---------------
def sech(x):
  return 2/(np.exp(x)+np.exp(-x))

def PTpot(xvalues):
  V0 = 0.15 # depth
  kappa = 0.1 # width
  return -V0 * sech(kappa*xvalues)**2

def zero_potential(xvalues):
    return np.zeros_like(xvalues)
# -------------------- now, do it ---------------
if __name__ == "__main__":
    endT = 400
    Nt = 800
    startX = -200
    endX = 200
    Nx = 800

    sigma = 5 # for gaussian pulse
    mu = -50
    ampl = 150
    width= 0.2 # for triangle pulse
    k = 1  # for square pulse

    deltat, timevalues, deltax, xvalues = gridmaker(endT,Nt,endX,Nx,startX)
    courant = c * deltat / deltax
    print("courant number = %.2f" % courant)

    ### potential
    potential = PTpot(xvalues)
    # potential = zero_potential(xvalues)
    plot_potential(xvalues,potential)

    Phi0, Pi0 = IVmaker('gauss',xvalues,sigma,mu,width,k,ampl) # phi: Zeilen: Zeit, Spalten: Ort
    bc = 'open_iii'
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential)
    # Phi, Pi = wave_evolution1D_4th_order(Phi0,Pi0,timevalues,xvalues,bc,potential)
    # Phi, Pi = wave_evolution1D_6th_order(Phi0,Pi0,timevalues,xvalues,bc,potential)
    print("calculation finished.")
    # Nt_plot = 7 # how many snap shots are plotted
    Etotal = total_energy(Phi,Pi)

    ### Plotting the results
    plot_energy_evolution(Etotal,timevalues)
    plot_xt_evolution_heatmap(timevalues,xvalues,Phi)
    xindex = 250
    plot_amplitude_evolution(timevalues,Phi[:,xindex],xvalues[xindex])
    plot_animation(xvalues, timevalues, Phi, Pi,'mp4')

    ### save as csv file
    # np.savetxt("results.csv", Phi, delimiter = ',', fmt = '%.6e')
