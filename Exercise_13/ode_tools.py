# -*- coding: utf-8 -*-
"""
ode_tools.py
"""

def euler_step(f,t,h,y,stages,get_metrics,payload_separates):
    '''
    Calculate one Euler step
    '''
    if get_metrics:
        state_dot,metrics = f(t,h,y,stages,get_metrics,payload_separates)
        return y + state_dot*h, metrics
    else:
        return y + f(t,y)*h

def midpoint_step(f,t,h,y,stages,get_metrics,payload_separates):
    '''
    Calculate one midpoint step
    '''
    y_mid = y + f(t,h,y,stages,get_metrics)*h/2
    
    if get_metrics:
        state_dot_mid,metrics = f(t+h/2,h,y_mid,stages,get_metrics,payload_separates)
        return y + state_dot_mid*h, metrics
    else:
        return y + f(t+h/2,y_mid)*h

def rk4_step(f,t,h,y,stages,get_metrics,payload_separates):
    '''
    Calculate one RK4 step
    '''
    if get_metrics:
        state_dot_1,metrics_1 = f(t,           h, y,                         stages, get_metrics, payload_separates)
        state_dot_2,metrics_2 = f(t + 0.5 * h, h, y + 0.5 * state_dot_1 * h, stages, get_metrics, payload_separates)
        state_dot_3,metrics_3 = f(t + 0.5 * h, h, y + 0.5 * state_dot_2 * h, stages, get_metrics, payload_separates)
        state_dot_4,metrics_4 = f(t + h,       h, y +       state_dot_3 * h, stages, get_metrics, payload_separates)
        
        state = y + h / 6.0 * (state_dot_1 + 2 * state_dot_2 + 2 * state_dot_3 + state_dot_4)
        metrics = (metrics_1 + 2 * metrics_2 + 2 * metrics_3 + metrics_4) / 6
        return state,metrics

    else:
        k1 = f(t,h,y, stages, get_metrics, payload_separates)
        k2 = f(t + 0.5 * h, h, y + 0.5 * k1 * h, stages, get_metrics, payload_separates)
        k3 = f(t + 0.5 * h, h, y + 0.5 * k2 * h, stages, get_metrics, payload_separates)
        k4 = f(t + h,       h, y +       k3 * h, stages, get_metrics, payload_separates)
    
        return y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

methods = {
    'euler': euler_step,
    'midpoint': midpoint_step,
    'rk4': rk4_step
}
