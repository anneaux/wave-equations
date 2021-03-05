import numpy as np

from finite_differences.example_functions import *

### constants
c = 1
coeffDict = {2:[1,-2,1],4:[-1/12,4/3,-5/2,4/3,-1/12], 6:[1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90],
8:[-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560]} # from https://doi.org/10.1090/S0025-5718-1988-0935077-0 bzw. https://en.wikipedia.org/wiki/Finite_difference_coefficient

### wave equation 0.2
# d phi / d t = phi
# d Pi / d t = c^2 d^2 phi / d x^2

### discretize time and space with uniform grid
def gridmaker1(endT,nt,endX,nx,startX=0,startT=0):
  dt = (endT-startT)/(nt)
  timevalues = np.linspace(startT,endT,nt+1)
  dx = (endX-startX)/(nx)
  xvalues = np.linspace(startX,endX,nx+1)
  print('dt = %.3f (%.f points), dx = %.3f (%.f points)' %(dt,nt,dx,nx))
  return dt,timevalues,dx,xvalues

def gridmaker(endT,maxX,courant=1):
  dt = 0.5
  dx = c * dt / courant

  # courant = c * deltat / deltax
  nt = int(endT/dt)
  timevalues = np.linspace(0,endT,nt+1)

  nx = int(2*maxX/dx)
  xvalues = np.linspace(-maxX,maxX,nx+1)
  print('dt = %.2f (%.f points), dx = %.2f (%.f points)' %(dt,nt,dx,nx))
  return dt,timevalues,dx,xvalues



def wave_evolution1D(phi0,Pi0,timevalues,xvalues,bc,potential,order=2):
  if bc != "periodic" and order != 2:
      print("for this accuracy order we only implemented 'periodic' boundary conditions, results may not be correct.")
  binomcoeffs = coeffDict[order]
  Nt = len(timevalues)
  Nx = len(xvalues)

  deltax = (xvalues[Nx-1]-xvalues[0])/(Nx-1)
  # print(deltax)
  deltat = timevalues[Nt-1]/(Nt-1)
  phi = np.zeros([Nt,Nx+order])     # add 1 ghost point at each end for each order
  Pi = np.zeros([Nt,Nx+order])
  ho = int(order/2) # half order, i.e. number of ghosts on each side

  ### first time step
  phi[0,ho:Nx+ho] = phi0          # fill inner points with given Phi0 and Pi0
  Pi[0,ho:Nx+ho] = Pi0
  potential = np.insert(potential,0,ho*[potential[ho]])   # expand potential to fit to ghosts
  potential = np.append(potential,ho*[potential[-1]])

  def rhs(phi,Pi,t):
    if bc == "periodic":
      phi[0:ho] = phi[-order-1:-ho-1] # filling ghost points at the beginning
      phi[-1:-ho-1:-1] = phi[order:order-ho:-1] # filling ghost points at the end
      Pi[0:ho] = Pi[-order-1:-ho-1]  # filling ghost points at the beginning
      Pi[-1:-ho-1:-1] = Pi[order:order-ho:-1] # filling ghost points at the end
    elif bc == "Dirichlet": # like a string
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
      phi[0] =  - deltax*(2*Pi[1] - Pi[2])/c + phi[1]
      phi[-1] =  - deltax*(2*Pi[-2] - Pi[-3])/c + phi[-2]
      Pi[0] =  (phi[0] - phi[1])/deltax #(2 * phi[1] - phi[2])/deltax + Pi[1] #
      Pi[-1] = (phi[-1] - phi[-2])/deltax # (2 * phi[-2] - phi[-3])/deltax + Pi[-2] # 

    # compute second spatial derivative (d^2 phi / dx^2) with FiniteDifferencing
    d2phidx2= np.zeros(Nx+order)
    for ix in range(ho,Nx+ho): # computing only inner points
      d2phidx2[ix] = 1/deltax**2 * sum(binomcoeffs[i]*phi[ix+i-ho] for i in range(0,order+1))
    dphidt = Pi
    dPidt = c**2 * d2phidx2 + potential*phi
    return dphidt, dPidt

  # for i in range(0,order+1):
  #     print("use index:", i-ho, "with coefficient ", binomcoeffs[i] )

  t = 1 #dummy value
  # time iteration (RK4 method)
  for i in range(0,Nt-1):
    k1_phi, k1_Pi = rhs(phi[i,:], Pi[i], t)
    k2_phi, k2_Pi = rhs(phi[i,:] + 0.5*deltat*k1_phi,Pi[i,:] + 0.5*deltat*k1_Pi,t + 0.5*deltat)
    k3_phi, k3_Pi = rhs(phi[i,:] + 0.5*deltat*k2_phi,Pi[i,:] + 0.5*deltat*k2_Pi ,t + 0.5*deltat)
    k4_phi, k4_Pi = rhs(phi[i,:] + deltat*k3_phi,Pi[i,:] + deltat*k3_Pi ,t + deltat)

    phi[i+1,:] = phi[i,:] + deltat*(1/6*k1_phi + 1/3*k2_phi +1/3*k3_phi + 1/6*k4_phi)
    Pi[i+1,:] = Pi[i,:] + deltat*(1/6*k1_Pi + 1/3*k2_Pi +1/3*k3_Pi + 1/6*k4_Pi)

  print("calculation finished.")
  return phi[:,ho:Nx+ho], Pi[:,ho:Nx+ho] # return only inner points

#--------------------- take a look at the energy -------------------
#-------------------------------------------------------------------
def energy(q,p):        #calculate energy from position q(phi) and inertia p(pi)
    m=1
    E = 0.5* p**2 / m + 0.5* q**2 *m
    return E
def total_energy(phi,pi,Nx):
    (rows,columns) = np.shape(phi)
    Etotal = np.zeros(rows)
    E = energy(phi,pi)
    for i in range(0,rows):   # for all times sum up individual energies
        #divide by number of columns to make E independent of Nx
        Etotal[i] = sum(E[i,1:Nx+1])/columns # do not consider ghost points
    Etotal = Etotal/max(Etotal)
    return Etotal


# -------------------- little helper function (InitialValuesMaker)-----------
def IVmaker(func,xvalues,sigma=1,mu=1,ampl=1,width=1,k=1):
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

def PT_potential(xvalues,depth,kappa):
  V0 = depth
  return -V0 * sech(kappa*xvalues)**2

def zero_potential(xvalues):
    return np.zeros_like(xvalues)
