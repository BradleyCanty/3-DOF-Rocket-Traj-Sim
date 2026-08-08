# -*- coding: utf-8 -*-
"""
ode_tools.py
"""

def euler_step(f,t,h,y,pitch_programs,get_metrics=False):
    '''
    Calculate one Euler step
    '''
    if get_metrics:
        state_dot,metrics = f(t,h,y,pitch_programs,get_metrics)
        return y + state_dot*h, metrics
    else:
        return y + f(t,y)*h

def midpoint_step(f,t,h,y,pitch_programs,get_metrics=False):
    '''
    Calculate one midpoint step
    '''
    y_mid = y + f(t,h,y,pitch_programs,get_metrics=False)*h/2
    
    if get_metrics:
        state_dot_mid,metrics = f(t + 0.5 * h, h, y_mid, pitch_programs,get_metrics=True)
        return y + state_dot_mid*h, metrics
    else:
        return y + f(t + 0.5 * h,y_mid)*h

def rk4_step(f,t,h,y,pitch_programs,get_metrics=False):
    '''
    Calculate one RK4 step
    '''
    if get_metrics:
        state_dot_1, metrics_1 = f(t,           h, y,                         pitch_programs, get_metrics=True)
        state_dot_2, metrics_2 = f(t + 0.5 * h, h, y + 0.5 * state_dot_1 * h, pitch_programs, get_metrics=True)
        state_dot_3, metrics_3 = f(t + 0.5 * h, h, y + 0.5 * state_dot_2 * h, pitch_programs, get_metrics=True)
        state_dot_4, metrics_4 = f(t + h,       h, y +       state_dot_3 * h, pitch_programs, get_metrics=True)
        
        state = y + h / 6.0 * (state_dot_1 + 2 * state_dot_2 + 2 * state_dot_3 + state_dot_4)
        metrics = (metrics_1 + 2 * metrics_2 + 2 * metrics_3 + metrics_4) / 6
        return state,metrics

    else:
        k1 = f(t,h,y)
        k2 = f(t + 0.5 * h, h, y + 0.5 * k1 * h)
        k3 = f(t + 0.5 * h, h, y + 0.5 * k2 * h)
        k4 = f(t + h,       h, y +       k3 * h)
    
        return y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

methods = {
    'euler': euler_step,
    'midpoint': midpoint_step,
    'rk4': rk4_step
}
