import numpy as np
# from rungekutta.rk4 import rk4
# from finite_differences.finite_differences_o2 import *
from finite_differences.example_functions import *
import matplotlib.pyplot as plt

### constants
c = 1

### wave equation 0.2
# d phi / d t = phi
# d Pi / d t = c^2 d^2 phi / d x^2

### discretize time and space with uniform grid
endT = 1
deltat = 0.01
Nt = int(endT/deltat)
timevalues = np.linspace(0,endT,Nt)

endX = 1
deltax = 0.01
Nx = int(endX/deltax)
xvalues = np.linspace(0,endX,Nx)
# print("Nx = %.1f" % Nx)

courant = c * deltat / deltax
print("courant number = %.2f" % courant)

def wave_evolution1D(phi0,Pi0,timevalues,xvalues):
  Nt = len(timevalues)
  Nx = len(xvalues)
  phi = np.zeros([Nt,Nx])
  Pi = np.zeros([Nt,Nx])

  phi[0,:] = phi0
  Pi[0,:] = Pi0

  def time_diff(phi, Pi, t): # where u = [phi, Pi]
    dphidt = Pi # because d/dt phi = Pi
    # d/dt Pi = c^2 * d2/dx2 phi, the latter shall be computed using FD (along the discretized x axis)
    d2phidx2= np.zeros(Nx)
    for i in range(1,Nx-1): # weil nur die inneren Werte genommen werden sollen, die Randwerte = 0
      d2phidx2[i] = 1/deltax**2 * (phi[i+1] - 2*phi[i] + phi[i-1])
    dPidt = c**2 * d2phidx2
    return dphidt, dPidt

  t = 1
  # time iteration (RK4 method)
  for i in range(0,Nt-1):
    k1_phi,k1_Pi  = time_diff(phi[i,:], Pi[i], t)
    k2_phi, k2_Pi = time_diff(phi[i,:] + 0.5*deltat*k1_phi,Pi[i,:] + 0.5*deltat*k1_Pi,t + 0.5*deltat)
    k3_phi,k3_Pi = time_diff(phi[i,:] + 0.5*deltat*k2_phi,Pi[i,:] + 0.5*deltat*k2_Pi ,t + 0.5*deltat)
    k4_phi,k4_Pi = time_diff(phi[i,:] + deltat*k3_phi,Pi[i,:] + deltat*k3_Pi ,t + deltat)

    phi[i+1,:] = phi[i,:] + deltat*(1/6*k1_phi + 1/3*k2_phi +1/3*k3_phi + 1/6*k4_phi)
    Pi[i+1,:] = Pi[i,:] + deltat*(1/6*k1_Pi + 1/3*k2_Pi +1/3*k3_Pi + 1/6*k4_Pi)

  return phi, Pi

# choose f_4, f_5 or g_a here (for the latter specify a = ...)
# phi0 = g_a(xvalues,0.025)
def gaussian(x,sigma,mu):
  return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(- (x-mu)**2/np.sqrt(2*sigma**2))
phi0 = gaussian(xvalues,0.005,0.5)

Pi0 = 0
phi, Pi = wave_evolution1D(phi0,Pi0,timevalues,xvalues)

############################################

# ### Plotting the solutions in a position-time diagram ###
# times = np.arange(0,Nt*deltat,deltat)
# plt.plot(times, phi_t)
# plt.xlabel('t')
# plt.ylabel('phi(t)')
# # # plt.axis('square')
# # # plt.xlim(-12,12)
# # # plt.ylim(49.94,50.008)
# # plt.legend()
# plt.grid(color = 'gainsboro')
# plt.savefig("plot-phi_t.pdf")

### Plotting the time evolution in a diagram with multiple lines
# from matplotlib import cm

# Nt_plot = 10 # how many snap shots are plotted
# Blues = cm.get_cmap('Blues_r',Nt_plot)

# for i in np.linspace(0,Nt-1,Nt_plot).astype(int):
#   plt.plot(xvalues, phi[i,:], label ="%.2f s" % timevalues[i], c = Blues(i/Nt))
# plt.xlabel('x')
# plt.ylabel('phi(x,t)')
# # # plt.axis('square')
# plt.xlim(0, max(xvalues))
# plt.ylim(np.amin(phi[:,:]),np.amax(phi[:,:]))
# plt.legend()
# plt.grid(color = 'gainsboro')
# plt.savefig("plot-phi(x,t).pdf")


### Plotting the time evolution in an animation ######################
import matplotlib
matplotlib.use('Agg')
import matplotlib.animation

def init_animation():
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

### write as gif
ani.save('WE-animation.gif', writer='imagemagick', fps=15)

### write as mp4
# Writer = matplotlib.animation.writers['ffmpeg'] # Set up formatting for the movie files
# mywriter = Writer(fps=15, metadata=dict(artist='AW'), bitrate=1800)
# ani.save('WE-animation.mp4', writer=mywriter)
