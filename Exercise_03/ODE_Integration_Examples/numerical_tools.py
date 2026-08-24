# -*- coding: utf-8 -*-
"""
numerical_tools.py
"""

import numpy as np
import ode_tools as ot

def propagate_ode(ode,state0,tspan,dt,method='rk4'):
    """
    This integrates the first order ordinary differential equations (specified
    by the 'ode' parameter) at each time step (specified by the 'tspan' and 'dt' 
    paramerters) using a numerical integration method (specified by the
    'method' parameter), and returns the times and associated state vectors

    Parameters
    ----------
    ode : function
        The first order ODE to be solved, i.e., the derivative of the state
    state0 : numpy array
        Specifies the initial state
    tspan : 1x2 array
        Contains the start time and the stop time
    dt : float
        the time step
    method : string
        Specifies the solver to use:
        'euler'    -> Eulers method
        'midpoint' -> midpoint method
        'rk4'      -> Runge Kutta 4 method

    Returns
    -------
    times : nx1 numpy array
        array of n times
    states : nxm numpy array
        array of m states, specified at n times
    """
    
    func        = ot.methods[method]
    times       = np.arange(tspan[0],tspan[1]+dt/2,dt)
    steps       = len(times)
    states      = np.zeros((steps, len(state0)))
    states[0]   = state0
    
    for step in range(steps - 1):
        states[step + 1] = func(ode,times[step], states[step], dt)
        
    return times, states




