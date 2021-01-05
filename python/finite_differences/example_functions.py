import numpy as np
from scipy import signal
############# define analytical functions ###########
def f_1(x):
    return (x-1/2)**2+x
def f_1_prime(x):
    return 2*x
def f_1_prime_prime(x):
    return 2+x*0
tex_f_1= r'$f_1=(x-\frac{1}{2})^{2}+x$'

def f_2(x):
    return (x-1/2)**3+(x-1/2)**2+x
def f_2_prime(x):
    return 3*(x**2)-x+(3/4)
def f_2_prime_prime(x):
    return 6*x-1
tex_f_2 = r'$f_2=(x-\frac{1}{2})^{3}+(x-\frac{1}{2})^{2}+x$'

def f_3(x):
    return np.sqrt(x)
def f_3_prime(x):
    return (1/2)*(x**(-1/2))
def f_3_prime_prime(x):
    return (-1/4)*(x**(-3/2))
tex_f_3 = r'$f_3=\sqrt{x}$'

def f_4(x,sigma = 1,mu=1,a=1):
    return np.sin(12*np.pi*x)
def f_4_prime(x,sigma = 1,mu=1,a=1):
    return np.cos(12*np.pi*x)*12*np.pi
def f_4_prime_prime(x,sigma = 1,mu=1,a=1):
    return -np.sin(12*np.pi*x)*(12*np.pi)**2
tex_f_4 = r'$f_4=\sin (12\pi x)$'

def f_5(x,sigma = 1,mu=1,a=1):
    return np.sin(12*np.pi*x)**4
def f_5_prime(x,sigma = 1,mu=1,a=1):
    return 4*12*np.pi*((np.sin(12*np.pi*x))**3)*(np.cos(12*np.pi*x))
def f_5_prime_prime(sigma = 1,mu=1,a=1):
    return 4*((12*np.pi)**2)*(3*(np.sin(12*np.pi*x)**2)*(np.cos(12*np.pi*x)**2)-np.sin(12*np.pi*x)**4)
tex_f_5 = r'$f_5=\sin^4 (12\pi x)$'

#set parameter for function g_a
def g_a(x,sigma = 1,mu=1,a=1):
    return np.exp(-a*x**2)
def g_a_prime(x,sigma = 1,mu=1,a=1):
    return np.exp(-a*x**2)*(-2*a*x)
def g_a_prime_prime(x,sigma = 1,mu=1,a=1):
    return np.exp(-a*x**2)*(4*a**2*x**2-2*a)
tex_g_a = r'$g_a=\exp(-ax^2)$'

### gaussian wave packet
def gaussian(x,sigma = 1,mu=1,a=1):
  return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(- (x-mu)**2/np.sqrt(2*sigma**2))
def gaussian_drv(x,sigma = 1,mu=1,a=1):
  return  -(x-mu)/(sigma**2 * np.sqrt(np.pi)) * np.exp(-(x-mu)**2/np.sqrt(2 * sigma**2))

### square pulse wave packet
def squares(x,k):
    return signal.square(2 * np.pi * k * (x-0.25))
def squares_drv(x,k):
    return np.zeros(len(x))
    
### triangle pulse wave packet
def f_triangle(xvalues,width,mu):
    T = np.zeros(len(xvalues))
    for i,x in enumerate(xvalues):
        T[i] = triangle(x,width,mu)
    return T
def f_triangle_drv(xvalues,width,mu):
    T = np.zeros(len(xvalues))
    for i,x in enumerate(xvalues):
        T[i] = triangle_drv(x,width,mu)
    return T
def triangle(x,width,mu):
    if x>mu: return triangle(1-x, width, mu)
    if x<mu-width:  return 0
    if x>mu-width:  return x-(mu-width)
    if x == mu: return 1
def triangle_drv(x,width,mu):
    if x>mu: return -1*triangle_drv(1-x, width, mu)
    if x<mu-width:  return 0
    if x>mu-width:  return 0.5
    if x == mu: return 0

func = [f_1,f_2,f_3,f_4,f_5,g_a]        #list of functions
func_1st_deriv = [f_1_prime,f_2_prime,f_3_prime,f_4_prime,f_5_prime,g_a_prime]      #list of 1st analyitcal function derivatives
func_2nd_deriv = [f_1_prime_prime,f_2_prime_prime,f_3_prime_prime,f_4_prime_prime,f_5_prime_prime,g_a_prime_prime] #list of 2nd analyitcal function derivatives

tex_f = [tex_f_1,tex_f_2,tex_f_3,tex_f_4,tex_f_5,tex_g_a]
