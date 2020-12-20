import numpy as np
from finite_differences.example_functions import *
import matplotlib.pyplot as plt
from scipy import signal

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
    elif bc == "open":
      # variant (i) WTF????
      # phi[0] = 2 * phi[1] - phi[2]
      # phi[-1] = 2 * phi[-2] - phi[-3]
      # Pi[0] = 2 * Pi[1] - Pi[2]
      # Pi[-1] = 2* Pi[-2] - Pi[-3]
      # variant (ii)
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


### gaussian wave packet
def gaussian(x,sigma,mu):
  return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(- (x-mu)**2/np.sqrt(2*sigma**2))
def gaussian_drv(x,sigma,mu):
  return  -(x-mu)/(sigma**2 * np.sqrt(np.pi)) * np.exp(-(x-mu)**2/np.sqrt(2 * sigma**2))

### square pulse wave packet
def squares(x,k):
    return signal.square(2 * np.pi * k * (x-0.25))
def squares_drv(x,k):
    return np.zeros(len(x))
### triangle pulse wave packet
def f_triangle(xvalues,width,mu):
    T = np.zeros(len(xvalues))
    for i,x in enumerate(xvalues):
        T[i] = triangle(x,width,mu)
    return T
def f_triangle_drv(xvalues,width,mu):
    T = np.zeros(len(xvalues))
    for i,x in enumerate(xvalues):
        T[i] = triangle_drv(x,width,mu)
    return T
def triangle(x,width,mu):
    if x>mu: return triangle(1-x, width, mu)
    if x<mu-width:  return 0
    if x>mu-width:  return x-(mu-width)
    if x == mu: return 1
def triangle_drv(x,width,mu):
    if x>mu: return -1*triangle_drv(1-x, width, mu)
    if x<mu-width:  return 0
    if x>mu-width:  return 0.5
    if x == mu: return 0
############################################

#### Plotting the solutions in a position-time diagram ###
def plot_xt(Nt, deltat, phi_t):
    times = np.arange(0,Nt*deltat,deltat)
    plt.plot(times, phi_t)
    plt.xlabel('t')
    plt.ylabel('phi(t)')
    # # plt.axis('square')
    # # plt.xlim(-12,12)
    # # plt.ylim(49.94,50.008)
    # plt.legend()
    plt.grid(color = 'gainsboro')
    plt.savefig("plots/plot-phi_t.png")

### Plotting the time evolution in a diagram with multiple lines
from matplotlib import cm
def plot_xt_evolution(timevalues,xvalues,phi,Nt_plot):
    Blues = cm.get_cmap('Blues_r',Nt_plot)

    for i in np.linspace(0,Nt-1,Nt_plot).astype(int):
      plt.plot(xvalues, phi[i,:], label ="%.2f s" % timevalues[i], c = Blues(i/Nt))
    plt.xlabel('x')
    plt.ylabel('phi(x,t)')
    # # plt.axis('square')
    plt.xlim(0, max(xvalues))
    plt.ylim(np.amin(phi[:,:]),np.amax(phi[:,:]))
    plt.legend()
    plt.grid(color = 'gainsboro')
    plt.savefig("plots/plot-phi(x,t).png")


### Plotting the time evolution in an animation ######################
import matplotlib
import matplotlib.animation
def plot_animation(xvalues, timevalues, phi, Pi):
  matplotlib.use('Agg')
  def init_animation(xvalues,phi):
    global line
    line, = ax.plot(xvalues, np.zeros_like(xvalues))
    ax.set_xlim(0, max(xvalues))
    # ax.set_ylim(0,10)
    ax.set_ylim(np.amin(phi[:,:]),np.amax(phi[:,:]))

  def animate(i):
    line.set_ydata(phi[i,:])
    timelabel.set_text('time: %.2f s' % timevalues[i])
    return line, timelabel

  fig3, ax3 = plt.subplots()
  ax3.set(xlabel = "x", ylabel = "phi(x)")
  line, = ax3.plot(xvalues, np.zeros_like(xvalues))
  # ax3.set_xlim(0, max(xvalues))
  ax3.set_ylim(np.amin(phi[:,:]),np.amax(phi[:,:]))
  # ax3.set_ylim(0,10)
  #, title = "time evolution of 1D wave")

  timelabel = ax3.text(0.02, 0.95, '', transform=ax3.transAxes)
  ani = matplotlib.animation.FuncAnimation(fig3, animate, frames=Nt, blit = True) #init_func=init_animation,

  # ### write as gif (gif stockt manchmal ein bisschen, geht außerdem sehr langsam zu speichern)
  # ani.save('plots/WE-animation.gif', writer='imagemagick', fps=15)

  ### write as mp4
  Writer = matplotlib.animation.writers['ffmpeg'] # Set up formatting for the movie files
  mywriter = Writer(fps=15, metadata=dict(artist='AW'), bitrate=1800)
  ani.save('plots/WE-animation.mp4', writer=mywriter)
#--------------------- take a look at the energy ------
def energy(q,p):        #calculate enervy from position q(phi) and inertia p(pi)
    m=1
    E = 0.5* p**2 / m + 0.5* q**2 *m
    return E
def total_energy(phi,pi):
    (rows,columns) = np.shape(phi)
    Etotal = np.zeros(rows)
    E = energy(phi,pi)
    for i in range(rows):
        #divide by number of columns to make E independent of Nx
        Etotal[i] = sum(E[i,:])/columns
    return Etotal
def plot_energy_evolution(Etotal,timevalues):
    fig, (ax1) = plt.subplots(1)
    ax1.plot(timevalues,Etotal, label='')
    ax1.set(xlabel='time $t$', ylabel='energy $E$')
    ax1.grid(color = 'gainsboro')
    plt.savefig('plots/WE_energy_evolution.png')

# -------------------- now, do it ---------------
if __name__ == "__main__":
    endT = 2
    Nt = 400
    endX = 1
    Nx = 100
    sigma = 0.005
    mu = 0.5
    width= 0.2
    k = 1

    deltat, timevalues, deltax, xvalues = gridmaker(endT,Nt,endX,Nx)
    # courant = c * deltat / deltax
    # print("courant number = %.2f" % courant)
# choose f_4, f_5, g_a (for latter specify a = ...) or gaussian here (for latter specify sigma and mu)

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
    Nt_plot = 7 # how many snap shots are plotted
    plot_energy_evolution(Etotal,timevalues)
    plot_xt_evolution(timevalues, xvalues, Phi, Nt_plot)
    plot_animation(xvalues, timevalues, Phi, Pi)
