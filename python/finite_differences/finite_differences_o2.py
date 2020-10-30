
################### first derivative ##############
def derivative_1_c(f,a,h=0.05):
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

    return (f(a + h) - f(a - h))/(2*h)
############
def derivative_2_c(f,a,h=0.05):
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

    return (f(a + h) - 2*f(a) + f(a - h))/(h**2)
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
