import matplotlib.pyplot as plt
import numpy as np
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
    Nxs = np.arange(Nx+1,(tests_h+1)*Nx+1,Nx)
    # print(np.shape(Nxs))
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
    Nxs = np.arange(Nx+1,(tests_cfl+1)*Nx+1,Nx)
    cfl = Nx/(Nt)
        # print(np.shape(Nxs))
    fig, (ax1) = plt.subplots(1)
    ax2 = ax1.twiny()
    ax1.plot(Nxs,abs_conv,'o-', label= 'absolute convergence' )
    ax1.plot(Nxs,self_conv,'o-', label= 'self convergence' )
    ax1.legend(title='CFL fixed: %.3f \n 1 period integration' %cfl)
    ax1.set(xlabel = 'number of spatial grid points', ylabel = '$Q= 2^p$' )
    ax1.grid(color = 'gainsboro')
    def tick_function(Nxs):
        Nts = (Nxs-1)/cfl +1
        return ["%d" % z for z in Nts]
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
