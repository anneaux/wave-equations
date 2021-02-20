# comparephi.py
from solver import *
from finite_differences.example_functions import *
from results_plotting import *
from waveconvergence_2 import convergence_test

def num_generator(cfl,Nt,bc,endT=1,endX=1):
    Nt *= endT
    Nx = int(cfl*Nt/endT)
    # num1
    dt,timevalues,dx,xvalues = gridmaker(endT, Nt, endX, Nx)
    potential = zero_potential(xvalues)
    Phi0, Pi0 = IVmaker('sine4',xvalues,1,1,1,1)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential)
    num1 = Phi[0::1]
    # print('timevalues:\n', timevalues )

    # print('timevalue for t=1', timevalues[int(Nt/endT)] )
    # print('num1(t= ',timevalues[int(Nt/endT)],' ):\n',num1[int(Nt/endT)])

    # num 2
    dt,timevalues,dx,xvalues = gridmaker(endT, 2*Nt, endX, 2*Nx)
    potential = zero_potential(xvalues)
    Phi0, Pi0 = IVmaker('sine4',xvalues,1,1,1,1)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential)
    num2 = Phi[0::2,0::2]
    # num3
    dt,timevalues,dx,xvalues = gridmaker(endT, 4*Nt, endX, 4*Nx)
    potential = zero_potential(xvalues)
    Phi0, Pi0 = IVmaker('sine4',xvalues,1,1,1,1)
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential)
    num3 = Phi[0::4,0::4]

    d2 = num1-num2
    d3 = num1-num3
    fig2 = plot_xt_evolution_heatmap(timevalues[0::4], xvalues[0::4], d2)
    fig3 = plot_xt_evolution_heatmap(timevalues[0::4], xvalues[0::4], d3)
    fig2.savefig('plots/convergence_2/differences/d2_endT_%s_nt_%s.png' %(endT,Nt))
    fig3.savefig('plots/convergence_2/differences/d3_endT_%s_nt_%s.png' %(endT,Nt))

    # ani2 = plot_animation(xvalues[0::4],timevalues[0::4],d2)
    # ani3 = plot_animation(xvalues[0::4],timevalues[0::4],d3)
    # Writer = matplotlib.animation.writers['ffmpeg'] # Set up formatting for the movie files
    # mywriter = Writer(fps=15, metadata=dict(artist='AW'), bitrate=1800)
    # ani2.save('plots/convergence_2/differences/d2_endT_%s_nt_%s.mp4' %(endT,Nt), writer=mywriter)
    # ani3.save('plots/convergence_2/differences/d3_endT_%s_nt_%s.mp4' %(endT,Nt), writer=mywriter)


if __name__ == "__main__":
    nt = 300
    cfl = 0.5
    endT = 2
    endX = 1

    bc = 'periodic'

    num_generator(cfl,nt,bc,endT,endX)
