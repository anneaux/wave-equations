import numpy as np
############# define analytical functions ###########
def f_1(x):
    return (x-1/2)**2+x
def f_1_prime(x):
    return 2*x
def f_1_prime_prime(x):
    return 2+x*0

def f_2(x):
    return (x-1/2)**3+(x-1/2)**2+x
def f_2_prime(x):
    return 3*x**2-x+(3/4)
def f_2_prime_prime(x):
    return 6*x-1

def f_3(x):
    return np.sqrt(x)
def f_3_prime(x):
    return (1/2)*x**(-1/2)
def f_3_prime_prime(x):
    return (-1/4)*x**(-3/2)

def f_4(x):
    return np.sin(12*np.pi*x)
def f_4_prime(x):
    return np.cos(12*np.pi*x)*12*np.pi
def f_4_prime_prime(x):
    return -np.sin(12*np.pi*x)*(12*np.pi)**2

def f_5(x):
    return np.sin(12*np.pi*x)**4
def f_5_prime(x):
    return 4*12*np.pi*((np.sin(12*np.pi*x))**3)*(np.cos(12*np.pi*x))
def f_5_prime_prime(x):
    return 4*((12*np.pi)**2)*(3*(np.sin(12*np.pi*x)**2)*(np.cos(12*np.pi*x)**2)-np.sin(12*np.pi*x)**4)

a = 5               #set parameter for function g_a
def g_a(x):
    return np.exp(-a*x**2)
def g_a_prime(x):
    return np.exp(-a*x**2)*(-2*a*x)
def g_a_prime_prime(x):
    return np.exp(-a*x**2)*(4*a**2*x**2-2*a)

func = [f_1,f_2,f_3,f_4,f_5,g_a]        #list of functions
func_1st_deriv = [f_1_prime,f_2_prime,f_3_prime,f_4_prime,f_5_prime,g_a_prime]      #list of 1st analyitcal function derivatives
func_2nd_deriv = [f_1_prime_prime,f_2_prime_prime,f_3_prime_prime,f_4_prime_prime,f_5_prime_prime,g_a_prime_prime] #list of 2nd analyitcal function derivatives
