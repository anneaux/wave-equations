import numpy as np

from finite_differences.example_functions import *
from results_plotting import *
### constants
c = 1
coeffDict = {2:[1,-2,1],4:[-1/12,4/3,-5/2,4/3,-1/12], 6:[1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90],
8:[-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560]} # from https://doi.org/10.1090/S0025-5718-1988-0935077-0 bzw. https://en.wikipedia.org/wiki/Finite_difference_coefficient

### wave equation 0.2
# d phi / d t = phi
# d Pi / d t = c^2 (d^2 phi / d x^2 + d^2 phi / d y^2)

### discretize time and space with uniform grid
def gridmaker(endT,nt,endX,nx,startX=0,startT=0):
  dt = (endT-startT)/(nt)
  timevalues = np.linspace(startT,endT,nt+1)
  dx = (endX-startX)/(nx)
  xvalues = np.linspace(startX,endX,nx+1)
  print('dt = %.2f (%.f points), dx = %.2f (%.f points)' %(dt,nt,dx,nx))
  return dt,timevalues,dx,xvalues


def gridmaker(endT,maxX,courant=1):
  dt = 0.5
  dx = c * dt / courant

  # courant = c * deltat / deltax
  nt = int(endT/dt)
  timevalues = np.linspace(0,endT,nt+1)

  nx = int(2*maxX/dx)
  xvalues = np.linspace(-maxX,maxX,nx+1)
  # print('dt = %.2f (%.f points), dx = %.2f (%.f points)' %(dt,nt,dx,nx))
  return dt,timevalues,dx,xvalues



def wave_evolution2D(phi0,Pi0,timevalues,xvalues,bc,potential,order):
  # if bc != "periodic" and order != 2:
  #     print("for this accuracy order we only implemented 'periodic' boundary conditions, results may not be correct.")
  binomcoeffs = coeffDict[order]
  Nt = len(timevalues)
  Nx = len(xvalues)
  Ny = len(xvalues)
  yvalues = xvalues

  deltax = (xvalues[Nx-1]-xvalues[0])/(Nx-1)
  # print(deltax)
  deltay = deltax
  deltat = timevalues[Nt-1]/(Nt-1)

  phi = np.zeros([Nt,Nx+order,Ny+order])  # add 1 ghost point at each end for each order
  Pi = np.zeros([Nt,Nx+order,Ny+order])
  ho = int(order/2) # half order, i.e. number of ghosts on each side

  ### first time step
  phi[0,ho:Nx+ho,ho:Ny+ho] = phi0          # fill inner points with given Phi0 and Pi0
  Pi[0,ho:Nx+ho,ho:Ny+ho] = Pi0

  ### expand initial distro to ghosts -> TODO this should be done better!
  phi[0,0,ho:Nx+ho] = phi0[0,:]
  phi[0,-1,ho:Nx+ho] = phi0[-1,:]
  phi[0,ho:Nx+ho,0] = phi0[:,0]
  phi[0,ho:Nx+ho,-1] = phi0[:,-1]

  Pi[0,0,ho:Nx+ho] = Pi0[0,:]
  Pi[0,-1,ho:Nx+ho] = Pi0[-1,:]
  Pi[0,ho:Nx+ho,0] = Pi0[:,0]
  Pi[0,ho:Nx+ho,-1] = Pi0[:,-1]


  # potential = np.insert(potential,0,ho*[potential[ho]])   # expand potential to fit to ghosts
  # potential = np.append(potential,ho*[potential[-1]])

  def rhs(phi,Pi,t):
    if bc == "open":         # variant (iii)
      ### bottom row
      ghostrow_b = np.zeros(Ny+order)
      for iy in range(ho,Ny+ho): #iterate over inner points only
        first = not (iy==ho)
        last = not (iy==Ny+ho-1)
        straightabove = 2*Pi[1,iy]-Pi[2,iy]
        diagleft = 2*Pi[1,iy-1]-Pi[2,iy-1]
        diagright = 2*Pi[1,iy+1]-Pi[2,iy+1]

        ghostrow_b[iy] = -deltax/c * (straightabove/2 + first*diagleft/4 + last*diagright/4) + phi[1,iy]/2 + last*phi[1,iy+1]/4 + first*phi[1,iy-1]/4

      ### top row
      ghostrow_t = np.zeros(Ny+order)
      for iy in range(ho,Ny+ho): #iterate over inner points only
        first = not (iy==ho)
        last = not (iy==Ny+ho-1)
        straightbelow = 2*Pi[-2,iy]-Pi[-3,iy]
        diagleft = 2*Pi[-2,iy-1]-Pi[-3,iy-1]
        diagright = 2*Pi[-2,iy+1]-Pi[-3,iy+1]

        ghostrow_t[iy] = -deltax/c * (straightbelow/2 + first*diagleft/4 + last*diagright/4) + phi[-2,iy]/2 + last*phi[-2,iy+1]/4 + first*phi[-2,iy-1]/4

      ### left column
      ghostcol_l = np.zeros(Nx+order)
      for ix in range(ho,Nx+ho): #iterate over inner points only
        first = not (ix==ho)
        last = not (ix==Nx+ho-1)
        right = 2*Pi[ix,1]-Pi[ix,2]
        diagabove = 2*Pi[ix+1,1]-Pi[ix+1,2]
        diagbelow = 2*Pi[ix-1,1]-Pi[ix-1,2]

        ghostcol_l[ix] = - deltay/c * (right/2 + last*diagabove/4 + first*diagbelow/4) + phi[ix,1]/2 + last*phi[ix+1,1]/4 + first*phi[ix-1,1]/4

      ### right column
      ghostcol_r = np.zeros(Nx+order)
      for ix in range(ho,Nx+ho): #iterate over inner points only
        first = not (ix==ho)
        last = not (ix==Nx+ho-1)
        left = 2*Pi[ix,-2]-Pi[ix,-3]
        diagabove = 2*Pi[ix+1,-2]-Pi[ix+1,-3]
        diagbelow = 2*Pi[ix-1,-2]-Pi[ix-1,-3]

        ghostcol_r[ix] = - deltay/c * (left/2 + last*diagabove/4 + first*diagbelow/4) + phi[ix,-2]/2 + last*phi[ix+1,-2]/4 + first*phi[ix-1,-2]/4

      phi[0,:] = ghostrow_b
      phi[-1,:] = ghostrow_t
      phi[:,0] = ghostcol_l
      phi[:,-1] = ghostcol_r


    # compute second spatial derivative (d^2 phi / dx^2) with FiniteDifferencing
    d2phidx2= np.zeros((Nx+order,Ny+order))
    d2phidy2= np.zeros((Nx+order,Ny+order))

    for ix in range(ho,Nx+ho): # computing only inner points
      d2phidx2[ix,:] = 1/deltax**2 * sum(binomcoeffs[i]*phi[ix+i-ho,:] for i in range(0,order+1))
      d2phidy2[:,ix] = 1/deltay**2 * sum(binomcoeffs[i]*phi[:,ix+i-ho] for i in range(0,order+1))

    dphidt = Pi
    dPidt = c**2 * (d2phidx2 + d2phidy2) #+ potential*phi
    return dphidt, dPidt

  # for i in range(0,order+1):
  #     print("use index:", i-ho, "with coefficient ", binomcoeffs[i] )

  ### time iteration (RK4 method)
  t = 1 #dummy value
  for i in range(0,Nt-1):
    k1_phi, k1_Pi = rhs(phi[i,:,:], Pi[i,:,:], t)
    k2_phi, k2_Pi = rhs(phi[i,:,:] + 0.5*deltat*k1_phi, Pi[i,:,:]+0.5*deltat*k1_Pi, t+0.5*deltat)
    k3_phi, k3_Pi = rhs(phi[i,:,:] + 0.5*deltat*k2_phi, Pi[i,:,:]+0.5*deltat*k2_Pi , t+0.5*deltat)
    k4_phi, k4_Pi = rhs(phi[i,:,:] + deltat*k3_phi, Pi[i,:,:]+deltat*k3_Pi , t+deltat)

    phi[i+1,:,:] = phi[i,:,:] + deltat*(1/6*k1_phi + 1/3*k2_phi +1/3*k3_phi + 1/6*k4_phi)
    Pi[i+1,:,:] = Pi[i,:,:] + deltat*(1/6*k1_Pi + 1/3*k2_Pi +1/3*k3_Pi + 1/6*k4_Pi)

  print("calculation finished.")
  return phi[:,ho:Nx+ho,ho:Ny+ho], Pi[:,ho:Nx+ho,ho:Ny+ho] # return only inner points



### gaussian wave packet
def gaussian(x,y,sigma=1,mux=0,muy=0,a=1):
# mu: mean value
# sigma: std deviation
  return a * np.exp(-(x-mux)**2/(2*sigma**2) - (y-muy)**2/(2*sigma**2)) # *1/np.sqrt(2*np.pi*sigma**2)
def gaussian_drv(x,y,sigma = 1,mu=1,a=1):
  return  a *(x-mu)/sigma**2 * gaussian(x,y,sigma,mu,a) # *1/np.sqrt(2*np.pi*sigma**2)

### plane wave front
def planewave(x,sigma=1,mux=0,a=1):
  # mu: mean value
  # sigma: std deviation
  return a * np.exp(-(x-mux)**2/(2*sigma**2))
def planewave_drv(x,sigma = 1, mux =0, a=1):
  return a *(x-mux)/sigma**2 * planewave(x,sigma,mux,a)




# ------------------- potential ---------------
def sech(x):
  return 2/(np.exp(x)+np.exp(-x))

def PT_potential(xvalues,yvalues,depth,kappa):
  V0 = depth
  return -V0 * sech(kappa*xvalues)**2


def zero_potential(xvalues):
    return np.zeros_like(xvalues)

# -------------------- now, do it ---------------
if __name__ == "__main__":

    sigma = 2 # for gaussian pulse
    mux = 0
    muy = 3
    ampl = 1

    depth = 0.15
    kappa = 0.1 # width

    ### numerical grid
    endT = 100
    maxX = 50
    courant = 1

    deltat, timevalues, deltax, xvalues = gridmaker(endT,maxX,courant)
    none, none, deltay, yvalues = gridmaker(endT,maxX,courant)
    Nx = len(xvalues)-1
    Ny = len(yvalues)-1
    Nt = len(timevalues)-1
    courant = c * deltat / deltax
    print('dt = %.2f (%.f points), dx = %.2f (%.f points)' %(deltat,Nt,deltax,Nx))
    print("courant number = %.2f" % courant)



    # -------------------- now, do it ---------------

    ### potential
    # potential = PT_potential(xvalues,yvalues, depth, kappa)
    potential = zero_potential(xvalues)
    # plot_potential(xvalues,potential)



### initial values
    Phi0 = np.zeros((Nx+1,Ny+1))
    Pi0 = np.zeros((Nx+1,Ny+1))

    ### plane wave
    # for ix in range(Nx+1):
    #   Phi0[ix,:] = planewave(xvalues[ix],sigma,mux,ampl)
    #   Pi0[ix,:] = planewave_drv(xvalues[ix],sigma,mux,ampl)

    ### gaussian blob
    for ix in range(Nx+1):
      for iy in range(Ny+1):
        Phi0[ix,iy] = gaussian(xvalues[ix],yvalues[iy],sigma,mux,muy,ampl)


    bc = 'open'
    order= 2
    Phi, Pi = wave_evolution2D(Phi0,Pi0,timevalues,xvalues,bc,potential,order)
    # print(Phi)

    plot_xt_evolution_heatmap(xvalues,yvalues,Phi[40,:,:])
    plot_2D_heatmap_animation(xvalues,yvalues,timevalues, Phi,'mp4')
