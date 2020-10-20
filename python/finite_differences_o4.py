################### first central derivative ##############
def derivative_1_c(f,a,h=0.05):
    '''Compute the difference formula for f'(a) with step size h.
    4th order accuracy
    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    h : number
        Step size in difference formula

    Returns
    -------
    float
        Difference formula:
            central: f(a+h) - f(a-h))/2h
    '''
    return (-f(a+2*h)+8*f(a + h) - 8*f(a - h)+ f(a-2*h) )/(12*h)

#################### second central derivative ###############
def derivative_2_c(f,a,h=0.05):
    '''Compute the difference formula for f''(a) with step size h.
    4th order accuracy
    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    h : number
        Step size in difference formula
    Returns
    -------
    float
        Difference formula:
            central: (f(a+h) - 2f(a) + f(a-h)) / h**2
    '''
    return (-f(a+2*h) + 16*f(a+h) - 30*f(a) + 16*f(a-h) - f(a-2*h)) / (12*(h**2))

################### first forward derivative ##############
def derivative_1_f(f,a,h=0.05):
    '''Compute the difference formula for f'(a) with step size h.
    4th order accuracy
    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    h : number
        Step size in difference formula

    Returns
    -------
    float
        Difference formula:
            central: (f(a+h) - f(a-h))/2h
    '''
    return (-25*f(a) + 48*f(a + h) -36*f(a + 2*h)+ 16*f(a+3*h) -3*f(a+4*h)) / (12*h)

#################### second forward derivative ###############
def derivative_2_f(f,a,h=0.05):
    '''Compute the difference formula for f''(a) with step size h.
    4th order accuracy
    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    h : number
        Step size in difference formula
    Returns
    -------
    float
        Difference formula:
            central: (f(a+h) - 2f(a) + f(a-h)) / h**2
    '''
    return ( 45*f(a) -154*f(a+h) +214*f(a+2*h) -156*f(a+3*h) +61*f(a+4*h) -10*f(a+5*h)) / (12*(h**2))

################### first backward derivative ##############
def derivative_1_b(f,a,h=0.05):
    '''Compute the difference formula for f'(a) with step size h.
    4th order accuracy
    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    h : number
        Step size in difference formula

    Returns
    -------
    float
        Difference formula:
            central: (f(a+h) - f(a-h))/2h
    '''
    return (25*f(a) - 48*f(a - h) +36*f(a - 2*h) - 16*f(a-3*h) +3*f(a-4*h)) / (12*h)

#################### second backward derivative ###############
def derivative_2_b(f,a, h=0.05):
    '''Compute the difference formula for f''(a) with step size h.
    4th order accuracy
    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    h : number
        Step size in difference formula
    Returns
    -------
    float
        Difference formula:
            central: (f(a+h) - 2f(a) + f(a-h)) / h**2
    '''
    return ( 45*f(a) -154*f(a-h) +214*f(a-2*h) -156*f(a-3*h) +61*f(a-4*h) -10*f(a-5*h)) / (12*(h**2))
