# -*- coding: utf-8 -*-
"""
ode_tools.py
"""

def euler_step(f,t,y,h):
    '''
    Calculate one Euler step
    '''
    return y + f(t,y)*h

def midpoint_step(f,t,y,h):
    '''
    Calculate one midpoint step
    '''
    y_mid = y + f(t,y)*h/2
    return y + f(t+h/2,y_mid)*h

def rk4_step(f,t,y,h):
    '''
    Calculate one RK4 step
    '''
    k1 = f(t,y)
    k2 = f(t + 0.5 * h, y + 0.5 * k1 * h)
    k3 = f(t + 0.5 * h, y + 0.5 * k2 * h)
    k4 = f(t + h,       y +       k3 * h)

    return y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

methods = {
    'euler': euler_step,
    'midpoint': midpoint_step,
    'rk4': rk4_step
}
