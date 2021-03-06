import numpy as np
import matplotlib.pyplot as plt

from finite_differences.example_functions import *

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
matplotlib.rcParams.update({'font.size': 20})

def plot_animation(xvalues, timevalues, phi, format = 'mp4'):
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
    ani.save('plots/WE_animation.mp4', writer=mywriter)
  print("...animation finished.")
  return ani

def plot_xt_evolution_heatmap(timevalues,xvalues,phi):
  from matplotlib.colors import LogNorm
  fig, ax = plt.subplots() #twilight_shifted
  im = ax.imshow(np.transpose(phi), cmap = 'twilight_shifted', origin = 'lower', extent = [min(timevalues),max(timevalues),min(xvalues),max(xvalues)])
  # norm=LogNorm(vmin=0.01, vmax=1)
  # maybe find better colormap? https://matplotlib.org/tutorials/colors/colormaps.html
  ax.set_aspect('auto')
  plt.xlabel('time')
  plt.ylabel('position')
  fig.colorbar(im)
  plt.savefig("plots/WE_phi_evolution_heatmap.png", bbox_inches = 'tight')
  # plt.show()
  return fig
def plot_xy_evolution_heatmap(timevalues,xvalues,phi):
  from matplotlib.colors import LogNorm
  fig, ax = plt.subplots()
  im = ax.imshow(np.transpose(phi), cmap = 'viridis', origin = 'lower', extent = [min(timevalues),max(timevalues),min(xvalues),max(xvalues)])
  # norm=LogNorm(vmin=0.01, vmax=1)
  # maybe find better colormap? https://matplotlib.org/tutorials/colors/colormaps.html
  ax.set_aspect('auto')
  plt.xlabel('x')
  plt.ylabel('y')
  fig.colorbar(im)
  plt.savefig("plots/WE_phi_evolution_heatmap.png", bbox_inches = 'tight')
  # plt.show()
  return fig
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
  plt.savefig('plots/PT_potential.png')

def plot_amplitude_evolution(timevalues,phi_at_xindex,x_at_xindex):
  fig, (ax1) = plt.subplots(1)
  ax1.plot(timevalues,phi_at_xindex, label='')
  ax1.set(xlabel='time', ylabel='phi')
  ax1.text(0.02, 0.95,"at x = %.2f" % x_at_xindex,transform=ax1.transAxes)
  # ax1.set_yscale('log')
  ax1.grid(color = 'gainsboro')
  plt.savefig('plots/WE_phi_evolution_onepoint.png')

def plot_amplitude_abs_evolution(timevalues,phi_at_xindex,x_at_xindex,maxphi):
  fig, (ax1) = plt.subplots(1)
  ax1.plot(timevalues,abs(phi_at_xindex), label='')
  ax1.set_ylim(0,maxphi)
  # ax1.set_yscale('log')
  ax1.set(xlabel='time', ylabel='abs(phi)')
  ax1.text(0.02, 0.95,"at x = %.2f" % x_at_xindex,transform=ax1.transAxes)
  ax1.grid(color = 'gainsboro')
  plt.savefig('plots/WE_phi_abs_evolution_onepoint_x%.1f.png' %(x_at_xindex))


def plot_amplitude_timestamp(xvalues,phi_at_tindex,t_at_tindex,depth,width):
  fig, (ax1) = plt.subplots(1)
  ax1.plot(xvalues,phi_at_tindex, label='')
  ax1.set(xlabel='position', ylabel='phi')
  ax1.text(0.02, 0.95,"at t = %.2f s" % t_at_tindex,transform=ax1.transAxes)
  # ax1.set_yscale('log')
  ax1.grid(color = 'gainsboro')
  # depth = 0.2
  # width = 0.5
  plt.savefig('plots/PT-timestamps-150s/WE_phi_timestamp_d%.2f_w%.2f.png' %(depth, width))



def plot_2D_heatmap_animation(xvalues,yvalues,timevalues, phi, format = 'mp4', filename = 'plots/WE_2D_animation.mp4'):
  print("start animation making...")
  matplotlib.use('Agg')
  Nt = len(timevalues)
  res_dpi=100
  figsize=4 #inch

   # animation function. This is called sequentially
  def animate(i):
      a = phi[i,:,:]
      # print(a)
      im.set_array(a)
      timelabel.set_text('time: %.2f s' % timevalues[i])
      return im, timelabel

  fig, ax = plt.subplots()
  ax.set(xlabel = "x", ylabel = "y")
  colmap = "nipy_spectral" # 
  # colmap = "gist_ncar"
  # colmap = "turbo"
  # colmap = "CMRmap"
  im = plt.imshow(phi[1,:,:],interpolation='none',cmap = colmap, origin = 'lower',extent = [min(xvalues),max(xvalues),min(yvalues),max(yvalues)])

  timelabel = ax.text(0.02, 0.95, '', transform=ax.transAxes)
  ani = matplotlib.animation.FuncAnimation(fig, animate, frames=Nt, blit = True)


  ### write as mp4
  if format == 'mp4':
    Writer = matplotlib.animation.writers['ffmpeg'] # Set up formatting for the movie files
    mywriter = Writer(fps=15, metadata=dict(artist='AW'), bitrate=1800)
    ani.save(filename, writer=mywriter)
  print("...animation finished.")
  return ani

def plot_2D_snapshot_heatmap(xvalues,yvalues,phi):
  from matplotlib.colors import LogNorm
  fig, ax = plt.subplots()
  im = ax.imshow(np.transpose(phi), cmap = 'viridis', origin = 'lower', extent = [min(xvalues),max(xvalues),min(yvalues),max(yvalues)])
  # norm=LogNorm(vmin=0.01, vmax=1)
  # maybe find better colormap? https://matplotlib.org/tutorials/colors/colormaps.html
  ax.set_aspect('auto')
  plt.xlabel('x')
  plt.ylabel('y')
  fig.colorbar(im)
  plt.savefig("plots/WE_2D_snapshot.png", bbox_inches = 'tight')
  # plt.show()
  return fig

from matplotlib.colors import Normalize
def plot_2D_snapshots_heatmaps(tvalues,xvalues,yvalues,phi,indices):
  # define circle 
  radius = 6.5
  theta = np.linspace(0,2*np.pi,100)
  circlex = radius * np.sin(theta)
  circley = radius *np.cos(theta)

  fig, axes = plt.subplots(1,len(indices)
    ,constrained_layout=True
    ,figsize=(20, 4.8))
  # cmap='viridis'
  cmap = "turbo"
  ext = [min(xvalues),max(xvalues),min(yvalues),max(yvalues)]
  normalizer=Normalize()#np.min(phi),np.max(phi))
  im=cm.ScalarMappable(norm=normalizer,cmap=cmap)
  for i,ax in enumerate(axes.flat):
    ax.imshow(np.transpose(phi[indices[i],:,:]),cmap=cmap,norm=normalizer,extent=ext)
    ax.set_title("t = %.1f" %(tvalues[indices[i]]))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.label_outer()
    ax.plot(circlex,circley,color ="gainsboro")
  fig.colorbar(im, ax=axes.ravel().tolist(),shrink=0.75)
  # plt.tight_layout()
  plt.savefig("plots/WE_2D_multi_snapshots.png", bbox_inches = 'tight')


def plot_scattered_waves(timevalues,phi_refl,x_refl,phi_trans,x_trans):
  fig, ax = plt.subplots(
        constrained_layout=True
        ,figsize=(6.4,4.8)
        )
  ax.plot(timevalues,abs(phi_refl), color = "tab:orange", label ="|φ| at x = -100")
  # ax.set_yscale('log')
  # ax.set_ylim(0,51)
  ax.set_xlabel('time')#, fontsize=fs)
  # ax.set_ylabel('|φ| at x = -100', color="tab:orange")#, fontsize=fs)
  # at x = %.2f" % x_refl
  ax.text(0.02, 0.9,"reflected wave",transform=ax.transAxes, color="tab:orange")
  ax.tick_params(axis="y",colors="tab:orange")
  # plt.yticks(fontsize=fs)
  # plt.xticks(fontsize=fs)
  ax.spines['left'].set_color('red')

  ax2=ax.twinx()
  ax2.plot(timevalues,abs(phi_trans), color = "tab:blue", label ="|φ| at x = 100")
  ax2.tick_params(axis="y",colors="tab:blue")
  # ax2.set_ylim(0,0.059)

  # ax2.set_ylabel('|φ| at x = 100', color="tab:blue")#, fontsize=fs)
  # at x = %.2f" % x_trans
  ax2.text(0.02, 0.8,"transmitted wave",transform=ax2.transAxes, color="tab:blue")

  ax.grid(color = 'gainsboro')
  # ax.legend()
  # ax2.legend()
  # plt.yticks(fontsize=fs)
  # plt.tight_layout()
  plt.savefig('plots/scattered_waves.png')

