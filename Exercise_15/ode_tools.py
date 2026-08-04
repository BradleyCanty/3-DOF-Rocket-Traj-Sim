# -*- coding: utf-8 -*-
"""
ode_tools.py
"""

def euler_step(ode,t,dt,state,stages,payload_separates):
    '''
    Calculate one Euler step
    '''
    state_dot,metrics = ode(t,dt,state,stages,payload_separates)
    return state + state_dot*dt, metrics

def midpoint_step(ode,t,dt,state,stages,payload_separates):
    '''
    Calculate one midpoint step
    '''
    state_mid = state + ode(t,dt,state,stages)*dt/2
    state_dot_mid,metrics = ode(t+dt/2,dt,state_mid,stages,payload_separates)
    return state + state_dot_mid*dt, metrics


def rk4_step(ode,t,dt,state,stages,payload_separates):
    '''
    Calculate one RK4 step
    '''
    state_dot_1,metrics_1 = ode(t,            dt, state,                          stages, payload_separates)
    state_dot_2,metrics_2 = ode(t + 0.5 * dt, dt, state + 0.5 * state_dot_1 * dt, stages, payload_separates)
    state_dot_3,metrics_3 = ode(t + 0.5 * dt, dt, state + 0.5 * state_dot_2 * dt, stages, payload_separates)
    state_dot_4,metrics_4 = ode(t + dt,       dt, state +       state_dot_3 * dt, stages, payload_separates)
    
    state = state + dt / 6.0 * (state_dot_1 + 2 * state_dot_2 + 2 * state_dot_3 + state_dot_4)
    metrics = (metrics_1 + 2 * metrics_2 + 2 * metrics_3 + metrics_4) / 6
    return state,metrics

methods = {
    'euler': euler_step,
    'midpoint': midpoint_step,
    'rk4': rk4_step
}
