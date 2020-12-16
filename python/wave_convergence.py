import numpy as np
from finite_differences.example_functions import *
import matplotlib.pyplot as plt
from solver import gridmaker, gaussian, gaussian_drv, wave_evolution1D
import math


def convergence(analytical, num_h, num_half_h):
    n = len(analytical)
    Q = [0]*n
    for i in range(n):
        Q[i] = abs(analytical[i]-num_h[i])/abs(analytical[i]-num_half_h[i])
    return Q

def self_convergence(num_h, num_half_h, num_quater_h):
    n = len(num_h)
    Q = [0]*n
    for i in range(n):
        Q[i] = abs(num_h[i]-num_half_h[i])/abs(num_half_h[i]-num_quater_h[i])
    return Q

def phi_select(phi,s):
    (a,b) = np.shape(phi) # a timepoint, b spacepoints
    original_length_phi = int(b/s)
    # print(original_length_phi)
    new_phi = np.zeros((a,int(b/s)))
    for i in range(int(b/s)):
        new_phi[:,i] = phi[:,s*i]
    return new_phi


def convergence_test(endT,Nt,endX,Nx,T,tests,startfunction,startfunction_prime):
    print('input: ',endT,Nt,endX,Nx,T,tests)
    # nt_conv:  Zeitpunkt in Zeitschritten, bei dem
    #           eine Periode vergangen ist
    nt_conv = int(Nt*T)
    # print('nt_conv = ', nt_conv)
    Phis = np.zeros((3,Nt,Nx))
# calculate 3 different phis vor 3 different
# spatial grid resolutions
    for k in range(0,3):
        # make grid
        # print('\n endT = %.1f, Nt = %d' %(endT,Nt))
        dt, timevalues, dx, xvalues = gridmaker(endT,Nt,endX,(2**k)*Nx)
        dt = endT/Nt
        # print('dt = %.3f , dx = %.3f' %(dt,dx))
        # print('shape of timevalues:', np.shape(timevalues))
        cour_N = c * dt / dx
        print("courant number = %.2f" % cour_N)
        # start conditions
        Phi0 = startfunction(xvalues)
        Pi0  = startfunction_prime(xvalues)
        # calculate
        Phi, Pi =  wave_evolution1D(Phi0,Pi0,timevalues,xvalues)
        Phi = phi_select(Phi,2**k)
        Phis[k,:,:] = Phi
# do the convergence tests
    # initialize arrays
    abs_conv = np.zeros((tests, Nx))
    self_conv = np.zeros((tests, Nx))
    for k in range(tests):
    # set stage for convergence tests
        analytical = f_4(np.linspace(0,endX,Nx))
        num_h = Phis[0, (k+1)*nt_conv, :]
        num_half_h = Phis[1, (k+1)*nt_conv, :]
        num_quater_h = Phis[2, (k+1)*nt_conv, :]
    # perform convergence tests
        abs_conv[k,:] = convergence(analytical, num_h, num_half_h)
        self_conv[k,:] = self_convergence(num_h, num_half_h, num_quater_h)
    return abs_conv, self_conv

def plot_convergence(abs_conv, self_conv, xvalues,tests):
    fig, (ax1,ax2) = plt.subplots(2,figsize=(10, 8))
    # fig, (ax1,ax2) = plt.subplots(2,figsize = (6,4))
    ax1.sharex(ax2)
    ax1.set_yscale('linear')
    ax2.set_yscale('linear')
    for i in range(tests):
        ax1.plot(xvalues,abs_conv[i,:],'-', label= '%d' %(i+1))
        ax2.plot(xvalues,self_conv[i,:],'-', label= '%d' %(i+1))
    ax1.set_title('absolute convergence')
    ax2.set_title('self convergence')#
    plt.xticks(np.arange(0, 1.1, 0.1))
    ax2.set_xlabel('x',loc='center')
    ax1.grid(color = 'gainsboro')
    ax2.grid(color = 'gainsboro')

    ax1.legend( title='periods integrated: ' )
    ax2.legend( title='periods integrated: ' )
    fig.savefig("wave-convergence.png")
    plt.show()



#--------------------------------------------------
#--------------------------------------------------

if __name__ == "__main__":
# set values for space time discretization
    c = 1
    endT = 1
    Nt = 6**3
    endX = 1
    Nx = 6**2
    T_per = 1/6 #T: Periodenlaenge
    n_tests = 5 # tests: Anzahl der Konvergenztests
    # test()

    dt, timevalues, dx, xvalues = gridmaker(endT,Nt,endX,Nx)

# choose f_4, f_5, g_a (for latter specify a = ...)
# or gaussian here (for latter specify sigma and mu)
    mu = 0.5
    sigma = mu/Nx
    a = 2
    Phi0 = f_5(xvalues)
    Pi0  = f_5_prime(xvalues)
    # phi0 = g_a(xvalues,a)#
    # Pi0 =  g_a_prime(xvalues,a)
    # phi0 = gaussian(xvalues,sigma,mu)
    # Pi0 = -c * gaussian_drv(xvalues,sigma,mu)

# run convergence tests
    abs_conv, self_conv = convergence_test(endT,Nt,endX,Nx,T_per,n_tests,f_4,f_4_prime)
# plot
    plot_convergence(abs_conv, self_conv, xvalues,n_tests)
    # print('Format Phi_0:', np.shape(Phi0))
    # print('Format Pi_0:', np.shape(Pi0))
    # print('Format Phi:', np.shape(Phi))
    # print('Format Pi:', np.shape(Phi))

# def test():
#     a = [1,2,3]
#     n1 = [1.1,1.5,2]
#     n2 = [1.01,1.8,2.5]
#     n3 = [1.011,1.9,2.6]
#     print(convergence(a,n1,n2), self_convergence(n1,n2,n3))
