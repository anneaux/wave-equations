import numpy as np
from finite_differences.example_functions import *
import matplotlib.pyplot as plt
from scipy import signal


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
    Nt = len(timevalues)
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

# my laptopresolution: 1920*1080
# 12,5" = 317.5mm, 16:9 ratio
# mylaptopDPI = 176
### Plotting the time evolution in an animation ######################
import matplotlib
import matplotlib.animation
def plot_animation(xvalues, timevalues, phi, Pi,format = 'mp4'):
  print("start animation making...")
  matplotlib.use('Agg')
  Nt = len(timevalues)
  res_dpi=100
  figsize=4 #inch
  # if len(xvalues) > figsize*res_dpi:
  #   phi = phi[:,0:-1:int(len(xvalues)/(res_dpi*figsize))]
  #   xvalues = xvalues[0:-1:int(len(xvalues)/(res_dpi*figsize))]
  def init_animation(xvalues,phi):
    global line
    line, = ax.plot(xvalues, np.zeros_like(xvalues))
    ax.set_xlim(0, max(xvalues))
    # ax.set_ylim(0,1000)
    # ax.set_yscale('log')

    # ax.set_ylim(np.amin(phi[:,:]),np.amax(phi[:,:]))

  def animate(i):
    line.set_ydata(phi[i,:])
    timelabel.set_text('time: %.2f s' % timevalues[i])
    return line, timelabel

  fig3, ax3 = plt.subplots()
  ax3.set(xlabel = "x", ylabel = "phi(x)")
  line, = ax3.plot(xvalues, np.zeros_like(xvalues)+0.01)
  # ax3.set_xlim(0, max(xvalues))
  ax3.set_ylim(np.amin(phi[:,:]),np.amax(phi[:,:]))
  # ax3.set_ylim(0,1000)
  # ax3.set_yscale('log')
  #, title = "time evolution of 1D wave")

  timelabel = ax3.text(0.02, 0.95, '', transform=ax3.transAxes)
  ani = matplotlib.animation.FuncAnimation(fig3, animate, frames=Nt, blit = True) #init_func=init_animation,

    ### write as gif (gif stockt manchmal ein bisschen, geht außerdem sehr langsam zu speichern)
  if format == 'gif':
    ani.save('plots/WE-animation.gif', writer='imagemagick', fps=15)

  ### write as mp4
  if format == 'mp4':
    Writer = matplotlib.animation.writers['ffmpeg'] # Set up formatting for the movie files
    mywriter = Writer(fps=15, metadata=dict(artist='AW'), bitrate=1800)
    ani.save('plots/WE-animation.mp4', writer=mywriter)
  print("...animation finished.")


def plot_xt_evolution_heatmap(timevalues,xvalues,phi):
  from matplotlib.colors import LogNorm
  fig, ax = plt.subplots()
  im = ax.imshow(np.transpose(phi), cmap = 'viridis', origin = 'lower', extent = [min(timevalues),max(timevalues),min(xvalues),max(xvalues)])
  # norm=LogNorm(vmin=0.01, vmax=1)
  # maybe find better colormap? https://matplotlib.org/tutorials/colors/colormaps.html
  ax.set_aspect('auto')
  plt.xlabel('time')
  plt.ylabel('position')
  plt.savefig("plots/WE_phi_evolution_heatmap.png", bbox_inches = 'tight')

def plot_energy_evolution(Etotal,timevalues):
    fig, (ax1) = plt.subplots(1)
    ax1.plot(timevalues,Etotal, label='')
    ax1.set(xlabel='time $t$', ylabel='energy $E$')
    ax1.grid(color = 'gainsboro')
    plt.savefig('plots/WE_energy_evolution.png')

def plot_potential(xvalues,potential):
  fig, (ax1) = plt.subplots(1)
  ax1.plot(xvalues,potential, label='')
  ax1.set(xlabel='x', ylabel='potential')
  ax1.grid(color = 'gainsboro')
  plt.savefig('plots/WE_PT_potential.png')

def plot_amplitude_evolution(timevalues,phi_at_x):
  fig, (ax1) = plt.subplots(1)
  ax1.plot(timevalues,phi_at_x, label='')
  ax1.set(xlabel='time', ylabel='phi')
  ax1.set_yscale('log')
  ax1.grid(color = 'gainsboro')
  plt.savefig('plots/WE_phi_evolution_onepoint.png')