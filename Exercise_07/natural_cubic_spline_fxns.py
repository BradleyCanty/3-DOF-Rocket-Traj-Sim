# -*- coding: utf-8 -*-
"""
natural_cubic_spline_fxns.py


Description:
    This constructs the cubic spline interpolant S for the function f, defined at the numbers
    x0 < x1 < ... < xn, satisfying S''(x0) = S''(xn) = 0
    
    INPUT:
    n; x0, x1, ... , xn; a0 = f(x0), a1 = f(x1), ... , an = f(xn)
    
    OUTPUT:
    aj, bj, cj, dj for j = 0, 1, ... , n - 1
    
    The coefficients are used in:
    S(x) = Sj(x) = aj + bj(x - xj) + cj(x - xj)^2 + dj(x-xj)^3 
    for xj <= x <= xj+1
    
Reference: Numerical Analysis by Burden and Faires, page 147

Example:
    Approximate the exponential f(x) = exp(x) using a natural spline and the following data points:
    (0,1), (1,e), (2,exp(2)), (3,exp(3))
"""
import math
import numpy as np

def get_cubic_spline_coeffs(x,y):
    '''

    Parameters
    ----------
    x : array of n independent values
    y : array of n dependent values

    Returns
    -------
    [a,b,c,d]: array of sub-arrays, where each sub-array has n-1 coefficients used in S(x) equation

    '''
    a = y
    h = [0] * (len(x)-1)
    for i in range(len(x)-1):
        h[i] = x[i+1] - x[i]
    
    alpha = [0] * (len(x) - 1)
    for i in range(1,len(x)-1):
        alpha[i] = 3/h[i]*(a[i+1] - a[i]) - 3/h[i-1]*(a[i]-a[i-1])
    
    l = [0] * (len(x))
    mu = [0] * (len(x))
    z = [0] * (len(x))
    l[0] = 1
    mu[0] = 0
    z[0] = 0 
    for i in range(1,len(x)-1):
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1] * mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i]
        
    b = [0] * len(x)
    c = [0] * len(x)
    d = [0] * len(x)
    
    for j in range(len(x)-2,-1,-1):
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    
    return [a[0:len(x)-1],b[0:len(x)-1],c[0:len(x)-1],d[0:len(x)-1]]

def get_cubic_spline_points(x,y,n):
    '''
    This fits a curve with n points between x and y
    '''
    #S(x) = Sj(x) = aj + bj(x - xj) + cj(x - xj)^2 + dj(x-xj)^3 
    #for xj <= x <= xj+1
    
    [a,b,c,d] = get_cubic_spline_coeffs(x,y)
    
    x_values = np.linspace(x[0],x[-1],n)
    y_values = np.zeros(len(x_values))
    
    for i in range(len(x_values)):
        for j in range(len(x)-1):
            if (x_values[i] >= x[j] and x_values[i] <= x[j+1]):
                y_values[i] = a[j] + b[j]*(x_values[i]-x[j]) + c[j]*(x_values[i]-x[j])**2 + d[j]*(x_values[i]-x[j])**3
    
    return x_values, y_values

def get_cubic_spline_point(x,y,x_value):
    '''
    This gets a single point using the curve fit to the x and y data
    '''
    #S(x) = Sj(x) = aj + bj(x - xj) + cj(x - xj)^2 + dj(x-xj)^3 
    #for xj <= x <= xj+1
    
    [a,b,c,d] = get_cubic_spline_coeffs(x,y)

    for j in range(len(x)-1):
        if (x_value >= x[j] and x_value < x[j+1]):
            y_value = a[j] + b[j]*(x_value-x[j]) + c[j]*(x_value-x[j])**2 + d[j]*(x_value-x[j])**3
            return y_value

if __name__ == "__main__":
    from matplotlib import pyplot as plt
    x = [0,1,2,3]
    y = [1,math.exp(1),math.exp(2),math.exp(3)]
    n = 101
    
    x_values,y_values = get_cubic_spline_points(x, y, n)

    y_values_actual = np.zeros(len(x_values))
    for i in range(len(x_values)):
        y_values_actual[i] = np.exp(x_values[i])
    
    plt.plot(x_values,y_values,'c',label='natural cubic spline fit')
    plt.plot(x,y,'b.')
    plt.plot(x_values,y_values_actual,'r--', label=r'$e^{x}$')
    plt.legend()
    plt.grid(True)

    #Get single value in the curve
    x_val = 1.75
    y_val = get_cubic_spline_point(x,y,x_val)

    #Plot the single point on the curve
    plt.plot(x_val,y_val,'y.',markersize=20)









