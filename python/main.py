from solver import *
from results_plotting import *


# -------------------- now, do it ---------------
if __name__ == "__main__":
    endT = 1
    Nt = 200
    startX = 0
    endX = 1
    Nx = 50

    sigma = 1/6 # for gaussian pulse
    mu = 1/2
    ampl = 1
    width= 0.2 # for triangle pulse
    k = 1  # for square pulse

    depth = 0.15
    kappa = 0.1 # width

# -------------------- now, do it ---------------
    deltat, timevalues, deltax, xvalues = gridmaker1(endT,Nt,endX,Nx,startX)
    courant = c * deltat / deltax
    print("courant number = %.2f" % courant)

    ### potential
    # potential = PT_potential(xvalues, depth, kappa)
    potential = zero_potential(xvalues)
    # plot_potential(xvalues,potential)

    Phi0, Pi0 = IVmaker('gauss',xvalues,sigma,mu,width,k,ampl) # phi: Zeilen: Zeit, Spalten: Ort
    bc = 'open_iii'
    order = 2
    Phi, Pi = wave_evolution1D(Phi0,Pi0,timevalues,xvalues,bc,potential,order)
    print("calculation finished.")

    # Etotal = total_energy(Phi,Pi,Nx)

    ### Plotting the results
    # plot_energy_evolution(Etotal,timevalues)
    plot_xt_evolution_heatmap(timevalues,xvalues,Phi)
    # xindex = 225
    # plot_amplitude_evolution(timevalues,Phi[:,xindex],xvalues[xindex])
    # tindex = 290
    # plot_amplitude_timestamp(xvalues,Phi[tindex,:],timevalues[tindex],depth,kappa)
    # plot_animation(xvalues, timevalues, Phi,'mp4')

    ### save as csv file
    # np.savetxt("results.csv", Phi, delimiter = ',', fmt = '%.6e')
