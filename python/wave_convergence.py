import numpy as np
from finite_differences.example_functions import *
import matplotlib.pyplot as plt
from solver import gridmaker, gaussian, gaussian_drv, wave_evolution1D
import math


def abs_convergence(analytical, num_h, num_half_h):
    print(np.shape(analytical))
    n = len(num_h)
    Q = 0
    Z = [0]*n
    N = [0]*n
    for i in range(n):
        Z[i] = (analytical[i]-num_h[i])**2
        N[i] = (analytical[i]-num_half_h[i])**2
    # print(np.sqrt(sum(Z)) , np.sqrt(sum(N)))
    Q = np.sqrt(sum(Z))/np.sqrt(sum(N))
    # print('abs conv ',Q)
    return Q

def self_convergence(num_h, num_half_h, num_quater_h):
    n = len(num_h)
    Q = 0
    Z = [0]*n
    N = [0]*n
    for i in range(n):
        Z[i] = abs(num_h[i]-num_half_h[i])**2
        N[i] = abs(num_half_h[i]-num_quater_h[i])**2
    Q = np.sqrt(sum(Z)) / np.sqrt(sum(N))
    return Q

def phi_select(phi,s):
    (a,b) = np.shape(phi) # a timepoint, b spacepoints
    original_length_phi = int(b/s)
    # print(original_length_phi)
    new_phi = np.zeros((a,int(b/s)))
    for i in range(int(b/s)):
        new_phi[:,i] = phi[:,s*i]
    return new_phi


def convergence_test_T_var(endT,Nt,endX,Nx,T,tests_T,startfunction,startfunction_prime):
# nt_conv:  Zeitpunkt in Zeitschritten, bei dem
#           eine Periode vergangen ist
    nt_conv = int(Nt*T)
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
        print("courant number = %.3f" % cour_N)
        # start conditions
        Phi0 = startfunction(xvalues)
        Pi0  = startfunction_prime(xvalues)
        # calculate
        Phi, Pi =  wave_evolution1D(Phi0,Pi0,timevalues,xvalues)
        Phi = phi_select(Phi,2**k)
        Phis[k,:,:] = Phi
# do the convergence tests
    # initialize
        abs_conv = np.zeros(tests_T)
        self_conv = np.zeros(tests_T)
    for k in range(tests_T):
    # select Phi at exact periods t = T
        analytical = startfunction(np.linspace(0,endX,Nx))
        num_h = Phis[0, (k+1)*nt_conv, :]
        num_half_h = Phis[1, (k+1)*nt_conv, :]
        num_quater_h = Phis[2, (k+1)*nt_conv, :]
    # perform convergence tests
        abs_conv[k] = abs_convergence(analytical, num_h, num_half_h)
        self_conv[k] = self_convergence(num_h, num_half_h, num_quater_h)
    return abs_conv, self_conv

def convergence_test_h_var(endT,Nt,endX,Nx,T_per,tests_h,startfunction,startfunction_prime):
    abs_conv = [0]*tests_h  # np.zeros((tests_h, Nx))
    self_conv = [0]*tests_h  # np.zeros((tests_h, Nx))
    for i in range(tests_h):
        abs_conv[i],self_conv[i] = convergence_test_T_var(endT,Nt,endX,(i+1)*Nx,T_per,1,startfunction,startfunction_prime)
    return abs_conv, self_conv

def convergence_test_cfl_var(endT,Nt,endX,Nx,T_per,tests_cfl,startfunction,startfunction_prime):
    abs_conv = [0]*tests_cfl  # np.zeros((tests_h, Nx))
    self_conv = [0]*tests_cfl  # np.zeros((tests_h, Nx))
    for i in range(tests_cfl):
        abs_conv[i],self_conv[i] = convergence_test_T_var(endT,(i+1)*Nt,endX,(i+1)*Nx,T_per,1,startfunction,startfunction_prime)
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

def plot_norms_T(abs_conv,self_conv,tests_T,Nx,Nt):
    Ts = range(1,tests_T+1)
    # print(np.shape(Ts))
    fig, (ax1) = plt.subplots(1)
    ax1.plot(Ts,abs_conv,'o-', label= 'absolute convergence' )
    ax1.plot(Ts,self_conv,'o-', label= 'self convergence' )
    ax1.legend(title='%d temporal grid points (CFL fixed) \n %d spatial grid points' %(Nt,Nx))
    ax1.set(xlabel = 'periods integrated', ylabel = '$Q= 2^p$')
    ax1.grid(color = 'gainsboro')
    plt.show()

def plot_norms_h(abs_conv,self_conv,tests_h,Nx,Nt):
    Nxs = np.arange(Nx,(tests_h+1)*Nx,Nx)
    print(Nxs)
    # print(np.shape(Ts))
    fig, (ax1) = plt.subplots(1)
    ax2 = ax1.twiny()
    ax1.plot(Nxs,abs_conv,'o-', label= 'absolute convergence' )
    ax1.plot(Nxs,self_conv,'o-', label= 'self convergence' )
    ax1.legend(title='%d temporal grid points (CFL varies) \n 1 period integration' %Nt)
    ax1.set(xlabel = 'number of spatial grid points', ylabel = '$Q= 2^p$')
    ax1.grid(color = 'gainsboro')
    def tick_function(Nxs):
        CFLS = Nxs/Nt
        return ["%.3f" % z for z in CFLS]
    ax1.set_xticks(Nxs)
    ax2.set_xticks(ax1.get_xticks())
    # ax2.set_xbound(ax1.get_xbound())
    ax2.set_xlim(ax1.get_xlim())
    ax1.tick_params(axis = 'y',direction = 'in')
    ax1.tick_params( axis = 'x', direction = 'in')
    ax2.tick_params(direction = 'in')
    ax2.set_xticklabels(tick_function(Nxs))
    ax2.set_xlabel('cfl number')
    plt.show()

def plot_norms_cfl(abs_conv,self_conv,tests_cfl,Nx,Nt):
    cfl = Nx/Nt
    Nxs = np.arange(Nx,(tests_cfl+1)*Nx,Nx)
    # print('plotting')
    # print('plot input: (abs_conv)', abs_conv, '(self_conv)', self_conv)
    # print('Nxs: ', Nxs)
    # print(np.shape(Nxs))
    fig, (ax1) = plt.subplots(1)
    ax2 = ax1.twiny()
    ax1.plot(Nxs,abs_conv,'o-', label= 'absolute convergence' )
    ax1.plot(Nxs,self_conv,'o-', label= 'self convergence' )
    ax1.legend(title='CFL fixed: %.3f \n 1 period integration' %cfl)
    ax1.set(xlabel = 'number of spatial grid points', ylabel = '$Q= 2^p$' )
    ax1.grid(color = 'gainsboro')
    def tick_function(Nxs):
        Nts = cfl * Nxs
        return ["%.3f" % z for z in Nts]
    ax1.set_xticks(Nxs)
    ax2.set_xticks(ax1.get_xticks())
    # ax2.set_xbound(ax1.get_xbound())
    ax2.set_xlim(ax1.get_xlim())
    ax1.tick_params(axis = 'y',direction = 'in')
    ax1.tick_params( axis = 'x', direction = 'in')
    ax2.tick_params(direction = 'in')
    ax2.set_xticklabels(tick_function(Nxs))
    ax2.set_xlabel('number of temporal grid points')
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
    n_tests_T = 5 # n_tests_T: Anzahl der perioden für die Konvergenz getestet wird
    n_tests_h = 10 # n_tests_h: Anzahl der x-Diskret. für die Konvergenz getestet wird
    n_tests_cfl = 10

    # dt, timevalues, dx, xvalues = gridmaker(endT,Nt,endX,Nx)
# choose f_4, f_5, g_a (for latter specify a = ...)
# or gaussian here (for latter specify sigma and mu)
    mu = 0.5
    sigma = mu/Nx
    a = 2
    # Phi0 = f_5(xvalues)
    # Pi0  = f_5_prime(xvalues)
    # phi0 = g_a(xvalues,a)#
    # Pi0 =  g_a_prime(xvalues,a)
    # phi0 = gaussian(xvalues,sigma,mu)
    # Pi0 = -c*gaussian_drv(xvalues,sigma,mu)
# run convergence test for grid resolution
    # abs_conv, self_conv = convergence_test_h_var(endT,Nt,endX,Nx,T_per,n_tests_h,f_4,f_4_prime)
    # plot_norms_h(abs_conv,self_conv,n_tests_h,Nx,Nt)
# run convergence test for up to 6 periods
    # abs_conv, self_conv = convergence_test_T_var(endT,Nt,endX,Nx,T_per,n_tests_T,f_4,f_4_prime)
    # plot_norms_T(abs_conv,self_conv,n_tests_T,Nx,Nt)
# run convergence test for CFL number
    abs_conv, self_conv = convergence_test_cfl_var(endT,Nt,endX,Nx,T_per,n_tests_cfl,f_4,f_4_prime)
    plot_norms_cfl(abs_conv,self_conv,n_tests_cfl,Nx,Nt)

    # print('Format Phi_0:', np.shape(Phi0))
    # print('Format Pi_0:', np.shape(Pi0))
    # print('Format Phi:', np.shape(Phi))
    # print('Format Pi:', np.shape(Phi))
