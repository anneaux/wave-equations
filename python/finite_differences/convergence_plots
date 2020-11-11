import numpy as np
import matplotlib.pyplot as plt

import finite_differences_o2 as dvo2
import finite_differences_o4 as dvo4
import example_functions as ex

################


func = ex.func                          #list of functions
func_1st_deriv = ex.func_1st_deriv      #list of 1st analyitcal function derivatives
func_2nd_deriv =  ex.func_2nd_deriv     #list of 2nd analyitcal function derivatives
tex_f = ex.tex_f                        #list of functions tex strings
############ set grid ###############################
n=20       #grid points
h=1/n       #grid spacing
x = np.linspace(0,1,n)

############ calculating #################
#function selection, input options = 0, 1, 2, 3, 4, 5(for g_a)
ex.a =1
fs = 5
#difference stencil method, input options= central, forward, backward
#m_1 = 'central'

y = func[fs](x)
y_prime = func_1st_deriv[fs](x)
y_prime_prime = func_2nd_deriv[fs](x)

dydx = dvo2.derivative_1_c(func[fs],x,h=h)    #calculate 1st derivative
dy2dx2 = dvo2.derivative_2_c(func[fs],x,h=h)    #calculate 2nd derivative
Delta_1 = y_prime - dydx              #analytical-difference 1st order
Delta_2 = y_prime_prime - dy2dx2        #analytical-difference 2nd order

################ plots ###################
#plt.rc('text', usetex=True)

fig, (ax1, ax2) = plt.subplots(2,figsize=(10, 9))

ax1.set_ylabel('')
#ax1.set_xlabel('x',loc='center')
ax1.plot(x,y,label='$f_%d(x)$' %(fs+1))
ax1.plot(x,y_prime,linestyle='dashed',label="Analytical 1st Derivative $f_%d^{(1)}(x)$" %(fs+1))
ax1.plot(x,dydx,marker = 'o',linestyle='dotted',label="Central 1st Difference $D^1_h f_%d(x)$" %(fs+1))
ax1.plot(x,y_prime_prime, linestyle='dashed',label="Analytical 2nd Derivative $f_%d^{(2)}(x)$" %(fs+1))
ax1.plot(x,dy2dx2,marker = 'o',linestyle='dotted',label="Central 2nd Difference $D^2_h f_%d(x)$" %(fs+1))
ax1.set_title('function and derivatives for ' '%s'  %(tex_f[fs]))
ax1.legend( title='grid size n=%d' % (n) )
ax1.grid(color = 'gainsboro')

ax2.set_ylabel('analytical - finite difference stencil')
ax2.set_xlabel('x',loc='center')
ax2.plot(x,Delta_1,marker = 'o', color = 'cyan', label=' 1st order')
ax2.plot(x,Delta_2, marker = 'o', color= 'hotpink', label=' 2nd order')
ax2.set_title('derivative difference')
ax2.legend()
ax2.grid(color = 'gainsboro')

plt.show()
