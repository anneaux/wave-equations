import numpy as np
import matplotlib.pyplot as plt

################### first derivative ##############
def derivative_1(f,a,method='central',h=0.05):
    '''Compute the difference formula for f'(a) with step size h.

    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    method : string
        Difference formula: 'forward', 'backward' or 'central'
    h : number
        Step size in difference formula

    Returns
    -------
    float
        Difference formula:
            central: f(a+h) - f(a-h))/2h
            forward: f(a+h) - f(a))/h
            backward: f(a) - f(a-h))/h
    '''
    if method == 'central':
        return (f(a + h) - f(a - h))/(2*h)
    elif method == 'forward':
        return (f(a + h) - f(a))/h
    elif method == 'backward':
        return (f(a) - f(a - h))/h
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.")

#################### second derivative ###############
def derivative_2(f,a,method='central',h=0.05):
    '''Compute the difference formula for f''(a) with step size h.

    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    method : string
        Difference formula: 'forward', 'backward' or 'central'
    h : number
        Step size in difference formula

    Returns
    -------
    float
        Difference formula:
            central: (f(a+h) - 2f(a) + f(a-h)) / h**2
            forward: (f(a+2h) - 2f(a+h) + f(a)) / h**2
            backward: (f(a) - 2f(a-h) + f(a-2h)) / h**2
    '''
    if method == 'central':
        return (f(a + h) - 2*f(a) + f(a - h))/(h**2)
    elif method == 'forward':
        return (f(a + 2*h) - 2*f(a + h) + f(a))/(h**2)
    elif method == 'backward':
        return (f(a) - 2*f(a - h) + f(a - 2*h))/(h**2)
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.")

############ set grid ###############################
n=40        #grid points
h=1/n       #grid spacing
x = np.linspace(0,1,n)
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

def f_a(x):
    return np.exp(-a*x**2)
def f_a_prime(x):
    return np.exp(-a*x**2)*(-2*a*x)
def f_a_prime_prime(x):
    return np.exp(-a*x**2)*(4*a**2*x**2-2*a)
a=1

############ calculating #################
y = f_a(x)
y_prime = f_a_prime(x)
y_prime_prime = f_a_prime_prime(x)
dydx = derivative_1(f_a,x,method = 'central',h=h)    #calculate 1st derivative
dy2dx2 = derivative_2(f_a,x,method = 'central',h=h)    #calculate 2nd derivative
Delta_1 = y_prime - dydx        #analytical-difference 1st order
Delta_2 = y_prime_prime - dy2dx2        #analytical-difference 2nd order

################ plots ###################
fig, (ax1, ax2) = plt.subplots(2,figsize=(10, 8))

ax1.set_ylabel('')
ax1.plot(x,y,label='y=f(x)')
ax1.plot(x,y_prime,linestyle='dashed',label="Analytical 1st Derivative y=f'(x)")
ax1.plot(x,dydx,linestyle='dotted',label="Central 1st Difference y=f'(x)")
ax1.plot(x,y_prime_prime,linestyle='dashed',label="Analytical 2nd Derivative y=f''(x)")
ax1.plot(x,dy2dx2,linestyle='dotted',label="Central 2nd Difference y=f''(x)")
ax1.legend()
ax1.grid(True)

ax2.set_ylabel('difference between derivates (analytical-difference)')
ax2.plot(x,Delta_1,color = 'cyan', label=' 1st order')
ax2.plot(x,Delta_2, color= 'hotpink', label=' 2nd order')
ax2.legend()
ax2.grid(True)

plt.show()
