from solver import *
from results_plotting import *

### incoming wave
sigma = 10 # for gaussian pulse
mu = -50
ampl = 50

### potential
depth = 0.15
kappa = 0.1 # width

### numerical grid
endT = 200
maxX = 150
courant = 1


# -------------------- now, do it ---------------
deltat, timevalues, deltax, xvalues = gridmaker(endT,maxX,courant)
print("courant number = %.2f" % courant)

potential = PT_potential(xvalues, depth, kappa)
# potential = zero_potential(xvalues)
# plot_potential(xvalues,potential)

Phi0, Pi0 = IVmaker('gauss',xvalues,sigma,mu,ampl) # phi: Zeilen: Zeit, Spalten: Ort

bc = 'open_iii'
order= 2
Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential,order)

Nx = len(xvalues)-1
Etotal = total_energy(Phi,Pi,Nx)

### Plotting the results
# plot_energy_evolution(Etotal,timevalues)
# plot_xt_evolution_heatmap(timevalues,xvalues,Phi)
# xindex = 100
# plot_amplitude_abs_evolution(timevalues,Phi[:,xindex],xvalues[xindex],47)
# xindex = 500
# plot_amplitude_abs_evolution(timevalues,Phi[:,xindex],xvalues[xindex],0.06)

plot_scattered_waves(timevalues,Phi[:,100],xvalues[100],Phi[:,500],xvalues[500])

# tindex = 290
# plot_amplitude_timestamp(xvalues,Phi[tindex,:],timevalues[tindex],depth,kappa)
# plot_animation(xvalues, timevalues, abs(Phi), 'mp4')
