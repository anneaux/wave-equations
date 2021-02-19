import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib.colors import LogNorm

from finite_differences.example_functions import *
from solver import *
from results_plotting import *
from convergence_plotting import *

def vonNeumann(phi,analytical):
    phi_err = analytical - phi
    (rows,columns) = np.shape(phi) #first index time, second space
    k = 2
    fouriers = np.zeros((rows,k*columns+1))
    for i in range(rows):
        fouriers[i,:] = fourier_coeff(phi_err[i,:],1)
        # fouriers[i,:] = np.fft.fft(phi_err[i,:],k*columns)
        # fouriers[i,:] = sp.fft.fft(phi_err[i,:])
        # print('fourier coefficients %d' %i,fouriers[i,:])
    G_factor = np.zeros((rows-1,k*columns+1))
    for i in range(0,rows-1):
        G_factor[i,:] = fouriers[i+1,:]/fouriers[i,:]
        # print('G_factor: %d' %i,G_factor[i,:])
    return G_factor

def fourier_coeff(y,period):
    c = np.zeros(2*len(y)+1)
    for n in range(-len(y),len(y)+1):
        c_m = 0
        for m in range(len(y)):
            c_m += y[m]*np.exp(-1j*2*np.pi*n*m/len(y))
        c[n] = c_m/c_m.size
    return c

#--------------------------
def vonNeumann_plot_heatmap(G_factor,endT,endX):

    fig, ax1 = plt.subplots()
    im = ax1.imshow(np.transpose(G_factor), cmap = 'viridis',
            origin = 'lower', extent=[0, endT, 0, endX], norm=LogNorm(vmin=0.01, vmax=10))
    fig.colorbar(im, ax=ax1)
    # maybe find better colormap? https://matplotlib.org/tutorials/colors/colormaps.html
    ax1.set(xlabel='time', ylabel='frequency')
    ### todo: set proper positions and times as ticks
    plt.savefig("plots/vonNeumann_heatmap.png", bbox_inches = 'tight')



########################

def test():
    endT = 1
    Nt = 50
    endX = 1
    Nx = 25
    # for gaussian pulse
    sigma = 0.005
    mu = 0.8
    ampl = 1
    # for triangle pulse
    width = 0.2
    # for square pulse
    k = 1

    deltat, timevalues, deltax, xvalues = gridmaker(endT,Nt,endX,Nx)
    # courant = c * deltat / deltax
    # print("courant number = %.2f" % courant)

    ### potential
    potential = np.zeros(np.shape(xvalues)) # PTpot(xvalues)
    # plot_potential(xvalues,potential)

    Phi0, Pi0 = IVmaker('sine4',xvalues,sigma,mu,width,k,ampl)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues, "periodic", potential)
    plot_xt_evolution_heatmap(timevalues,xvalues,Phi)

    ana = np.zeros(np.shape(Phi))
    for i in range(Nx):
        for j in range(Nt):
            ana[j,i] = sine_12_pi(i,j)
    G_factor = vonNeumann(Phi,ana)
    G_factor[0,:] = 0

    vonNeumann_plot_heatmap(G_factor,endT,endX)

test()
