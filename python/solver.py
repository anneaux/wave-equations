import numpy as np
from finite_differences.example_functions import *
import matplotlib.pyplot as plt

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


def wave_evolution1D(phi0,Pi0,timevalues,xvalues):
  Nt = len(timevalues)
  Nx = len(xvalues)
  deltax = xvalues[Nx-1]/Nx
  deltat = timevalues[Nt-1]/Nt
  phi = np.zeros([Nt,Nx+2])
  Pi = np.zeros([Nt,Nx+2])
  ### first time step
  phi[0,1:Nx+1] = phi0
  ### fill ghost points with values at boundary
  phi[0,0] = phi0[0]
  phi[0,Nx+1] = phi0[Nx-1]

  Pi[0,1:Nx+1] = Pi0
  ### fill ghost points with values at boundary
  Pi[0,0] = Pi0[0]
  Pi[0,Nx+1] = Pi0[Nx-1]

  def time_diff(phi, Pi, t): # where u = [phi, Pi]
    dphidt = Pi # because d/dt phi = Pi
    # d/dt Pi = c^2 * d2/dx2 phi, the latter shall be computed using FD (along the discretized x axis)
    d2phidx2= np.zeros(Nx+2)
    for ix in range(1,Nx+1): # computing only inner points
      d2phidx2[ix] = 1/deltax**2 * (phi[ix+1] - 2*phi[ix] + phi[ix-1])
    d2phidx2[-1] = d2phidx2[1]
    d2phidx2[0] = d2phidx2[-2]
    dPidt = c**2 * d2phidx2
    return dphidt, dPidt

  t = 1
  # time iteration (RK4 method)
  for i in range(0,Nt-1):
    # right boundary
    phi[i, -1] = phi[i,1]
    Pi[i, -1] = Pi[i,1]
    # left boundary
    phi[i, 0] = phi[i,-2]
    Pi[i, 0] = Pi[i,-2]


    k1_phi, k1_Pi  = time_diff(phi[i,:], Pi[i], t)
    k2_phi, k2_Pi = time_diff(phi[i,:] + 0.5*deltat*k1_phi,Pi[i,:] + 0.5*deltat*k1_Pi,t + 0.5*deltat)
    k3_phi, k3_Pi = time_diff(phi[i,:] + 0.5*deltat*k2_phi,Pi[i,:] + 0.5*deltat*k2_Pi ,t + 0.5*deltat)
    k4_phi, k4_Pi = time_diff(phi[i,:] + deltat*k3_phi,Pi[i,:] + deltat*k3_Pi ,t + deltat)

    phi[i+1,:] = phi[i,:] + deltat*(1/6*k1_phi + 1/3*k2_phi +1/3*k3_phi + 1/6*k4_phi)
    Pi[i+1,:] = Pi[i,:] + deltat*(1/6*k1_Pi + 1/3*k2_Pi +1/3*k3_Pi + 1/6*k4_Pi)

  return phi[:,1:Nx+1], Pi[:,1:Nx+1] # return only inner points


### gaussian wave packet
def gaussian(x,sigma,mu):
  return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(- (x-mu)**2/np.sqrt(2*sigma**2))
def gaussian_drv(x,sigma,mu):
  return  -(x-mu)/(sigma**2 * np.sqrt(np.pi)) * np.exp(-(x-mu)**2/np.sqrt(2 * sigma**2))

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
    plt.savefig("plot-phi_t.png")

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
    plt.savefig("plot-phi(x,t).png")


### Plotting the time evolution in an animation ######################
import matplotlib
import matplotlib.animation
def plot_animation(xvalues, timevalues, phi, Pi):
    matplotlib.use('Agg')
    def init_animation(xvalues,phi):
      global line
      line, = ax.plot(xvalues, np.zeros_like(xvalues))
      ax.set_xlim(0, max(xvalues))
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
    #, title = "time evolution of 1D wave")

    timelabel = ax3.text(0.02, 0.95, '', transform=ax3.transAxes)
    ani = matplotlib.animation.FuncAnimation(fig3, animate, frames=Nt, blit = True) #init_func=init_animation,

    ### write as gif (gif stockt manchmal ein bisschen, geht außerdem sehr langsam zu speichern)
    ani.save('WE-animation.gif', writer='imagemagick', fps=15)

    ### write as mp4
    Writer = matplotlib.animation.writers['ffmpeg'] # Set up formatting for the movie files
    mywriter = Writer(fps=15, metadata=dict(artist='AW'), bitrate=1800)
    ani.save('WE-animation.mp4', writer=mywriter)

# -------------------- now, do it ---------------
if __name__ == "__main__":
    endT = 1
    Nt = 100
    endX = 1
    Nx = 100
    deltat, timevalues, deltax, xvalues = gridmaker(endT,Nt,endX,Nx)
    courant = c * deltat / deltax
    # print("courant number = %.2f" % courant)
# choose f_4, f_5, g_a (for latter specify a = ...) or gaussian here (for latter specify sigma and mu)
    phi0 = f_4(xvalues)
    Pi0  = f_4_prime(xvalues)
    # phi0 = g_a(xvalues,20)#
    # Pi0 = 3*np.zeros(len(phi0))#- g_a_prime(xvalues,20)
    # phi0 = gaussian(xvalues,0.005,0.5)
    # Pi0 = -c * gaussian_drv(xvalues,0.005,0.5)
    phi, Pi = wave_evolution1D(phi0,Pi0,timevalues,xvalues)

    Nt_plot = 5 # how many snap shots are plotted
    plot_xt_evolution(timevalues,xvalues,phi,Nt_plot)
    # plot_animation(xvalues, timevalues, phi, Pi)
