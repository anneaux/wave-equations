import numpy as np
from finite_differences.example_functions import *
import matplotlib.pyplot as plt
from solver import *
import math
from results_plotting import *
from convergence_plotting import *


def abs_convergence(analytical, num_h, num_half_h):
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
    (a,b) = np.shape(phi)       # a: timepoints, b: spacepoints
    orig_length_phi = int((b-1)/s+1)
    # print(orig_length_phi)
    new_phi = np.zeros((a,orig_length_phi))
    for i in range(orig_length_phi):
        new_phi[:,i] = phi[:,s*i]
    return new_phi


def convergence_test_T_var(endT,Nt,endX,Nx,T,tests_T,func,func_prime,sigma,mu,a,bc):
# nt_conv:  Zeitpunkt in Zeitschritten, bei dem
#           eine Periode vergangen ist
    nt_conv = int(Nt*T)
    Phis = np.zeros((3,Nt,Nx+1))
# calculate 3 different phis vor 3 different
# spatial grid resolutions
    for k in range(0,3):
        # make grid
        # print('\n endT = %.1f, Nt = %d' %(endT,Nt))
        dt, timevalues, dx, xvalues = gridmaker(endT,Nt,endX,(2**k)*Nx+1)
        dt = endT/Nt
        # print('dt = %.3f , dx = %.3f' %(dt,dx))
        # print('shape of timevalues:', np.shape(timevalues))
        # print('shape of xvalues:', np.shape(xvalues))
        cour_N = c * dt / dx
        print("courant number = %.3f" % cour_N)
        # start conditions
        Phi0 = func(xvalues,sigma,mu,a)
        Pi0  = func_prime(xvalues,sigma,mu,a)
        potential = zero_potential(xvalues)
        # calculate
        Phi, Pi =  wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential)
        Phi = phi_select(Phi,2**k)
        Phis[k,:,:] = Phi
# do the convergence tests
    # initialize
        abs_conv = np.zeros(tests_T)
        self_conv = np.zeros(tests_T)
    for k in range(tests_T):
    # select Phi at exact periods t = k*T
        analytical = func(np.linspace(0,endX,Nx+1))
        num_h = Phis[0, (k+1)*nt_conv, :]
        num_half_h = Phis[1, (k+1)*nt_conv, :]
        num_quater_h = Phis[2, (k+1)*nt_conv, :]
    # perform convergence tests
        print(np.len(xvalues),np.len(timevalues))
        abs_conv[k] = abs_convergence(analytical, num_h, num_half_h)
        self_conv[k] = self_convergence(num_h, num_half_h, num_quater_h)
    return abs_conv, self_conv

def convergence_test_h_var(endT,Nt,endX,Nx,T_per,tests_h,func,func_prime,sigma,mu,a,bc):
    abs_conv = [0]*tests_h  # np.zeros((tests_h, Nx))
    self_conv = [0]*tests_h  # np.zeros((tests_h, Nx))
    for i in range(tests_h):
        abs_conv[i],self_conv[i] = convergence_test_T_var(endT,Nt,endX,(i+1)*Nx,T_per,1,
                            func,func_prime,sigma,mu,a,bc)
    return abs_conv, self_conv

def convergence_test_cfl_var(endT,Nt,endX,Nx,T_per,tests_cfl,func,func_prime,sigma,mu,a,bc):
    abs_conv = [0]*tests_cfl  # np.zeros((tests_h, Nx))
    self_conv = [0]*tests_cfl  # np.zeros((tests_h, Nx))
    for i in range(tests_cfl):
        abs_conv[i],self_conv[i] = convergence_test_T_var(endT,(i+1)*Nt,endX,(i+1)*Nx,T_per,1,
                            func,func_prime,sigma,mu,a,bc)
    return abs_conv, self_conv


#--------------------------------------------------
#--------------------------------------------------

if __name__ == "__main__":
# set values for space time discretization
    c = 1
    z = 1   # helper vairable for convergence over integrated periods
    endT = z
    Nt = 2*z*6**3
    endX = 1
    Nx = 2*6**2
    T_per = 1/(z*6) #T: Periodenlaenge
    n_tests_T = z*6 - 1 # n_tests_T: Anzahl der perioden für die Konvergenz getestet wird
    n_tests_h = 10 # n_tests_h: Anzahl der x-Diskretisierungen bei konstanter CFL für die Konvergenz getestet wird
    n_tests_cfl = 10

    mu = 0.5
    sigma = mu/Nx
    a = 2
# run convergence test for up to 6 periods
    # abs_conv, self_conv = convergence_test_T_var(endT,Nt,endX,Nx,T_per,n_tests_T,
    #                     f_4,f_4_prime,sigma,mu,a,'open_iii')
    # plot_norms_T(abs_conv,self_conv,n_tests_T,Nx,Nt)
# run convergence test for increasing grid resolution (cfl increasing!)
    # abs_conv, self_conv = convergence_test_h_var(endT,Nt,endX,Nx,T_per,n_tests_h,
    #                     f_4,f_4_prime,sigma,mu,a,'periodic')
    # plot_norms_h(abs_conv,self_conv,n_tests_h,Nx,Nt)
# run convergence test for increasing grid resolution (cfl fixed!!!)
    abs_conv, self_conv = convergence_test_cfl_var(endT,Nt,endX,Nx,T_per,n_tests_cfl,
                        f_4,f_4_prime,sigma,mu,a,'open_iii')
    plot_norms_cfl(abs_conv,self_conv,n_tests_cfl,Nx,Nt)
