import numpy as np
import matplotlib.pyplot as plt

import finite_differences_o2 as dvo2
import finite_differences_o4 as dvo4
import example_functions as ex

################


func = [ex.f_1,ex.f_2,ex.f_3,ex.f_4,ex.f_5,ex.g_a]        #list of functions
func_1st_deriv = [ex.f_1_prime,ex.f_2_prime,ex.f_3_prime,ex.f_4_prime,ex.f_5_prime,ex.g_a_prime]      #list of 1st analyitcal function derivatives
func_2nd_deriv = [ex.f_1_prime_prime, ex.f_2_prime_prime, ex.f_3_prime_prime,ex.f_4_prime_prime,ex.f_5_prime_prime,ex.g_a_prime_prime] #list of 2nd analyitcal function derivatives

############ set grid ###############################
n=60       #grid points
h=1/n       #grid spacing
x = np.linspace(0,1,n)

############ calculating #################
fs = 3               #function selection, input options = 0, 1, 2, 3, 4, 5(for g_a)
method = 'central'

y = func[fs](x)
y_prime = func_1st_deriv[fs](x)
y_prime_prime = func_2nd_deriv[fs](x)

dydx = dvo4.derivative_1_c(func[fs],x,h=h)    #calculate 1st derivative
dy2dx2 = dvo4.derivative_2_c(func[fs],x,h=h)    #calculate 2nd derivative
Delta_1 = y_prime - dydx              #analytical-difference 1st order
Delta_2 = y_prime_prime - dy2dx2        #analytical-difference 2nd order

################ plots ###################
fig, (ax1, ax2) = plt.subplots(2,figsize=(10, 8))

ax1.set_ylabel('')
#ax1.set_xlabel('x',loc='center')
ax1.plot(x,y,label='y=f(x)')
ax1.plot(x,y_prime,linestyle='dashed',label="Analytical 1st Derivative y=f'(x)")
ax1.plot(x,dydx,marker = 'o',linestyle='dotted',label="Central 1st Difference y=f'(x)")
ax1.plot(x,y_prime_prime, linestyle='dashed',label="Analytical 2nd Derivative y=f''(x)")
ax1.plot(x,dy2dx2,marker = 'o',linestyle='dotted',label="Central 2nd Difference y=f''(x)")
ax1.set_title('function and derivatives for f%d' %(fs+1))
ax1.legend( title='grid size n=%d, %s stencil' % (n,method) )
ax1.grid(color = 'gainsboro')

ax2.set_ylabel('analytical - finite difference stencil')
ax2.set_xlabel('x',loc='center')
ax2.plot(x,Delta_1,marker = 'o', color = 'cyan', label=' 1st order')
ax2.plot(x,Delta_2, marker = 'o', color= 'hotpink', label=' 2nd order')
ax2.set_title('derivative difference')
ax2.legend()
ax2.grid(color = 'gainsboro')

plt.show()
