# -*- coding: utf-8 -*-
"""
math_fxns_ex2.py
"""
import numpy as np

def get_rocket_state(ry0, vy0, tt, soln_method):
    """
    Establish the ODEs to solve, then solve them numerically
    
    Parameters:
    -----------
    y0 : float array
        initial y position
    vy0 : float array
        2-tuple of initial velocities
    tt : float array
        array of times at which to evaluate the ODEs
    soln_method : string
        solution method for solving system of equations
    
    Returns:
    --------
    array of floats
        y positions, y velocities
    """
    
    r_earth = 6371000 #m, mean earth radius
    m_empty = 4000 #kg, empty mass
    m_propellant = 8800 #kg
    m_dot = 129.4 #kg/s, mass flow rate of fuel
    v_exhaust = 2050 #m/s, exhaust velocity
    t_burnout = m_propellant/m_dot #seconds
    g = 9.81
    
    def sys_of_ODEs(t,ry,vy):
        '''
        Establish system of ODEs to solve
        
        Parameters:
        -----------
        t: double
            time at which to evaluate the function
        y: double
            y position
        vy: double
            y velocity
        
        Returns:
        -------
        float array
            y velocity and y acceleration      
        '''
        
        if (t < t_burnout):
            m_inst = m_empty + m_propellant - m_dot*t #kg, instantaneous mass
            ay = v_exhaust * m_dot / m_inst - g
        else:
            m_inst = m_empty
            ay = -g
        return vy, ay, m_inst
    
    #Numerically integrate the system of ODEs
    if soln_method == 'Euler': #using tangential method, a.k.a. "Euler method" or "1st order Runge-Kutta method"
        ry = np.ones(len(tt)) * r_earth
        vy = np.zeros(len(tt))
        ay = np.zeros(len(tt))
        m = np.zeros(len(tt))
        
        ry[0] = ry0
        vy[0] = vy0
        
        dt = tt[1] - tt[0] #assume time step in supplied time array is constant
        
        for ti in range(0,len(tt)-1):
            drydt, dvydt, m_inst = sys_of_ODEs(tt[ti],ry[ti],vy[ti])
            ay[ti] = dvydt
            m[ti] = m_inst
            vy[ti+1] = dvydt*dt + vy[ti]
            ry[ti+1] = drydt*dt + ry[ti]
            if ry[ti+1] < r_earth:
                tt = tt[0:ti+1]
                ry = ry[0:ti+1]
                vy = vy[0:ti+1]
                ay = ay[0:ti+1]
                m = m[0:ti+1]
                break
        
        return tt, ry, vy, ay, m
    
     

    
    