import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

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
        fouriers[i,:] = fourier_coeff(phi[i,:],1)
        # fouriers[i,:] = np.fft.fft(phi_err[i,:],k*columns)
        # fouriers[i,:] = sp.fft.fft(phi_err[i,:])
        # print('fourier coefficients %d' %i,fouriers[i,:])
    G_factor = np.zeros((rows-1,k*columns+1))
    for i in range(0,rows-1):
        G_factor[i,:] = fouriers[i+1,:]/fouriers[i,:]
        print('G_factor: %d' %i,G_factor[i,:])
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
def vonNeumann_plot_heatmap(G_factor):
    fig, ax1 = plt.subplots()
    im = ax1.imshow(np.transpose(G_factor), cmap = 'viridis', origin = 'lower')
    fig.colorbar(im, ax=ax1)
    # maybe find better colormap? https://matplotlib.org/tutorials/colors/colormaps.html
    plt.xlabel('time')
    plt.ylabel('position')
    ### todo: set proper positions and times as ticks
    plt.savefig("plots/vonNeumann_heatmap.png", bbox_inches = 'tight')

def vonNeumann_plot_two():
    return 0

########################

def test():
    endT = 1
    Nt = 50
    endX = 1
    Nx = 25
    # for gaussian pulse
    sigma = 0.005
    mu = 0.8
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

    Phi0, Pi0 = IVmaker('sine',xvalues,sigma,mu,width,k)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues, "open_iii", potential)

    ana = np.zeros(np.shape(Phi))
    for i in range(Nx):
        for j in range(Nt):
            ana[j,i] = sine_12_pi(i,j)
    G_factor = vonNeumann(Phi,ana)
    G_factor[0,:] = 0

    vonNeumann_plot_heatmap(G_factor)

test()
