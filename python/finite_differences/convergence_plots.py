import numpy as np
import matplotlib.pyplot as plt

import finite_differences_o2 as dvo2
import finite_differences_o4 as dvo4
import example_functions as ex

################
#list of functions
func = ex.func
#list of 1st analyitcal function derivatives
func_1st_deriv = ex.func_1st_deriv
#list of 2nd analyitcal function derivatives
func_2nd_deriv =  ex.func_2nd_deriv
#list of functions tex strings
tex_f = ex.tex_f
################## function selection #########################
#input options = 0, 1, 2, 3, 4, 5(for g_a)
ex.a = 1
fs = 2
########### set grid ###############################
n = 40      #grid points
x = np.linspace(0,1,n)

############ calculating #################
y = func[fs](x)
y_prime = func_1st_deriv[fs](x)
y_prime_prime = func_2nd_deriv[fs](x)
t=3
A = np.zeros((t,n))
h = [0]*3
dydx = A
dy2dx2 = A
Delta_prime = A
Delta_prime_prime = A
r = A

for a in range(t):
    #grid spacing for derivatives
    h[a]=(1/n)*(1/2)**a
    #calculate 1st derivative
    dydx[a] = dvo2.derivative_1_c(func[fs],x,h=h[a])
    #calculate 2nd derivative
    dy2dx2[a] = dvo2.derivative_2_c(func[fs],x,h=h[a])
################### outer convergence test ################
for a in range(t):
    r[a] = y_prime - dydx[a]

Q = abs(r[0])/(4*abs(r[1]))
################### self convergence test #################
T = (dydx[0]-dydx[1])/((dydx[1]-dydx[2]))
############################# plots #######################
#plt.rc('text', usetex=True)

fig, ax = plt.subplots(figsize=(10, 5))

ax.set_ylabel('')
ax.set_xlabel('x',loc='center')
ax.plot(x,T, linestyle= 'dotted',
    marker = 'o',
    label='$[D^1_h-D^1_{h/2}] / [D^1_{h/2}-D^1_{h/4}] $')
ax.set_title('self-convergence for ' '%s'  %(tex_f[fs]))
ax.legend( title='grid size n=%d' % (n) )
ax.grid(color = 'gainsboro')

plt.show()
