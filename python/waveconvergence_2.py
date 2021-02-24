# waveconvergence_2.py
from solver import *
from finite_differences.example_functions import *
from results_plotting import *
#------------------- plotting -------------

def self_convergence_plotting_2(cfl, Nt, selfQ, endT):
    fig, (ax1) = plt.subplots(1)
    for i in range(len(Nt)):
        timevalues = np.linspace(0,endT,int(Nt[i])+1)
        ax1.plot(timevalues,selfQ[i],'-',
                label = '%d timepoints' %( int(Nt[i])+1 ) )
    plt.xticks(np.linspace(0,endT,7))
    plt.title('self convergence')
    ax1.legend(title = 'cfl=%.2f' %(cfl))
    ax1.set(xlabel = 'time', ylabel = '$p=log_2(Q)$')
    ax1.grid(color = 'gainsboro')
    # plt.show()
    return fig
def abs_convergence_plotting_2(cfl, Nt, absQ, endT):
    fig, (ax1) = plt.subplots(1)
    for i in range(len(Nt)):
        timevalues = np.linspace(0,endT,int(Nt[i])+1)
        ax1.plot(timevalues,absQ[i],'-',
                label = '%d timepoints' %( int(Nt[i])/endT +1 ) )
    plt.xticks(np.linspace(0,endT,endT*6+1))
    plt.title('absolute convergence')
    ax1.legend(title = 'cfl=%.2f' %(cfl))
    ax1.set(xlabel = 'time', ylabel = '$p=log_2(Q)$')
    ax1.grid(color = 'gainsboro')
    # plt.show()
    return fig

#----------------------- conv tests -----------

def self_convergence(num1,num2,num3):
    p = np.zeros(np.shape(num1)[0])
    Q = np.zeros(np.shape(num1)[0])
    # print('shape self Q: ', np.shape(Q))
    for i in range(np.shape(num1)[0]):
        # difference1 = np.sqrt(np.sum(abs(num1[i,:]**2 - num2[i,:]**2)))
        # difference2 = np.sqrt(np.sum(abs(num2[i,:]**2 - num3[i,:]**2)))
        difference1 = np.sqrt(np.sum( (num1[i,:] - num2[i,:])**2 ))
        difference2 = np.sqrt(np.sum( (num2[i,:] - num3[i,:])**2 ))
        Q[i] = difference1/difference2
        p[i] = np.log2(difference1/difference2)
        # print('diff1: \n', difference1)
    # print('Q: \n',Q)
    return p
def abs_convergence(ana,num1,num2):
    p = np.zeros(np.shape(num1)[0])
    Q = np.zeros(np.shape(num1)[0])
    for i in range(np.shape(num1)[0]):
        # difference1 = np.sqrt(np.sum(abs(init**2 - num1[i,:]**2)))
        # difference2 = np.sqrt(np.sum(abs(init**2 - num2[i,:]**2)))
        difference1 = np.sqrt(np.sum( (ana - num1[i,:])**2 ))
        difference2 = np.sqrt(np.sum( (ana - num2[i,:])**2 ))
        Q[i] = difference1/difference2
        p[i] = np.log2(difference1/difference2)
    return p

#---------------main function -------------------------------------------

def convergence_test(cfl,Nt,bc,funchandle,endT=1,endX=1):
    Nx = int(cfl*Nt/endT)

    # num1
    dt,timevalues,dx,xvalues = gridmaker1(endT, Nt, endX, Nx)
    # print('xvalues\n',xvalues)
    potential = zero_potential(xvalues)
    Phi0, Pi0 = IVmaker(funchandle,xvalues,1,1,1,1)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential)
    num1 = Phi[0::1]
    # print('num1(t=0):\n',num1[0])
    # print('timevalues:\n', timevalues )
    # print('timevalue for t=1', timevalues[int(Nt/endT)] )
    # print('num1(t= ',timevalues[int(Nt/endT)],' ):\n',num1[int(Nt/endT)])
    # print('shape(num1): ', np.shape(num1))
    # # print('shape(timevalues1): ',np.shape(timevalues[0::1]))
    # print('timevalues1:', timevalues)

    # ana
    ana = np.zeros((Nt+1,Nx+1))
    for i in range(Nt+1):
        for j in range(Nx+1):
            ana[i,j] = sine_12_pi(xvalues[j],timevalues[i])

    # num 2
    dt,timevalues,dx,xvalues = gridmaker1(endT, 2*Nt, endX, 2*Nx)
    potential = zero_potential(xvalues)
    Phi0, Pi0 = IVmaker(funchandle,xvalues,1,1,1,1)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential,)
    num2 = Phi[0::2,0::2]
    # print('shape(num2): ', np.shape(num2))
    # # print('shape(timevalues2): ',np.shape(timevalues[0::2]))
    # print('timevalues2:\n', timevalues)

    # num3
    dt,timevalues,dx,xvalues = gridmaker1(endT, 4*Nt, endX, 4*Nx)
    potential = zero_potential(xvalues)
    Phi0, Pi0 = IVmaker(funchandle,xvalues,1,1,1,1)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential)
    num3 = Phi[0::4,0::4]
    # print('shape(num2): ', np.shape(num2))
    # print('timevalues3:', timevalues)

    selfQ = self_convergence(num1, num2, num3)
    absQ = abs_convergence(ana,num1,num2)
    return selfQ, absQ

if __name__ == "__main__":
    cfl = 0.5
    endT = 1
    endX = 1
    tests = 1
    bc = 'open_iii'
    Nt_step = 100
    k = 0         # k+Nt_step: starting value for Nts
    Nt = np.linspace(endT*(k+Nt_step),endT*(k+Nt_step*tests),tests)
    selfQs = []
    absQs = []
    fhandle = "sine4"
    for i in range(len(Nt)):
        selfQ, absQ = convergence_test(cfl,int(Nt[i]),bc,fhandle,endT,endX)
        selfQs.append(selfQ)
        absQs.append(absQ)
    fig1 = self_convergence_plotting_2(cfl, Nt, selfQs, endT)
    fig2 = abs_convergence_plotting_2(cfl, Nt, absQs, endT)
    fig1.savefig('plots/convergence_2/self_conv_%s.png' %(bc))
    fig2.savefig('plots/convergence_2/abs_conv_%s.png' %(bc))
