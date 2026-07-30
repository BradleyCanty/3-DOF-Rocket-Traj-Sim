# -*- coding: utf-8 -*-
"""
math_fxns_ex4.py
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
        solution method for solving system of equations: either Euler, Midpoint, or RK4
    
    Returns:
    --------
    array of floats
        y positions, y velocities
    """
    
    r_earth = 6371000 #m, mean earth radius
    
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
        m_empty = 4000 #kg, empty mass
        m_propellant = 8800 #kg
        m_dot = 129.4 #kg/s, mass flow rate of fuel
        v_exhaust = 2050 #m/s, exhaust velocity
        t_burnout = m_propellant/m_dot #seconds
        G = 6.6743E-11 #m^3 kg^-1 s^-2 
        m_earth = 5.97219E24 #kg
        g = lambda ry: G * m_earth / ry**2
        rho_air = 1.22 #kg/m^3
        D = 1.65 #m, diameter of V2 rocket
        A = np.pi*(D/2)**2
        Cd = .125 #drag coefficient of V2 rocket
        
        if (t < t_burnout):
            ay = v_exhaust * m_dot / (m_empty + m_propellant - m_dot*t) - g(ry) - (1/2 * rho_air * vy * abs(vy) * A * Cd) / (m_empty + m_propellant - m_dot*t)
        else:
            ay = -g(ry) - (1/2 * rho_air * vy * abs(vy) * A * Cd) / m_empty
        return vy, ay
    
    #Numerically integrate the system of ODEs
    if soln_method == 'Euler': #using tangential method, a.k.a. "Euler method" or "1st order Runge-Kutta method"
        ry = np.ones(len(tt)) * r_earth
        vy = np.zeros(len(tt))
        ay = np.zeros(len(tt))
        
        ry[0] = ry0
        vy[0] = vy0
        
        dt = tt[1] - tt[0] #assume time step in supplied time array is constant
        
        for ti in range(0,len(tt)-1):
            drydt, dvydt = sys_of_ODEs(tt[ti],ry[ti],vy[ti])
            ay[ti] = dvydt
            vy[ti+1] = dvydt*dt + vy[ti]
            ry[ti+1] = drydt*dt + ry[ti]
            if ry[ti+1] < r_earth:
                tt = tt[0:ti+1]
                ry = ry[0:ti+1]
                vy = vy[0:ti+1]
                ay = ay[0:ti+1]
                break
        
        return tt, ry, vy, ay
    
    elif soln_method == 'midpoint': #using midpoint method, a.k.a. "2nd order Runge-Kutta method"
        ry = np.ones(len(tt)) * r_earth
        vy = np.zeros(len(tt))
        
        ry[0] = ry0
        vy[0] = vy0
        
        dt = tt[1] - tt[0] #assume time step in supplied time array is constant
    
        for i in range(0,len(tt)-1):
            drydt, dvydt = sys_of_ODEs(tt,(ry[i],vy[i]))
            ryMid = drydt*dt/2 + ry[i]
            vyMid = dvydt*dt/2 + vy[i]
            
            drxdtMid, drydtMid, dvxdtMid, dvydtMid = sys_of_ODEs(tt,(ryMid,vyMid))
            ry[i+1] = drydtMid*dt + ry[i]
            vy[i+1] = dvydtMid*dt + vy[i]
        
        return tt, ry, vy
    
    elif soln_method == 'RK4': #using 4th order Runge-Kutta method
        ry = np.zeros(len(tt))
        vy = np.zeros(len(tt))
        
        ry[0] = ry0
        vy[0] = vy0
        
        dt = tt[1] - tt[0] #assume time step in supplied time array is constant
        
        for i in range(0,len(tt)-1):
            drydt1, dvydt1 = sys_of_ODEs(tt,ry[i],vy[i])
            ry2 = drydt1 * dt/2 + ry[i]
            vy2 = dvydt1 * dt/2 + vy[i]
            drydt2, dvydt2 = sys_of_ODEs(tt,ry2,vy2)
            ry3 = drydt2 * dt/2 + ry[i]
            vy3 = dvydt2 * dt/2 + vy[i]
            drydt3, dvydt3 = sys_of_ODEs(tt,ry3,vy3)
            ry4 = drydt3 * dt + ry[i]
            vy4 = dvydt3 * dt + vy[i]
            drxdt4, drydt4, dvxdt4, dvydt4 = sys_of_ODEs(tt,ry4,vy4)
            
            ry[i+1] = ry[i] + dt*(drydt1 + 2*drydt2 + 2*drydt3 + drydt4)/6
            vy[i+1] = vy[i] + dt*(dvydt1 + 2*dvydt2 + 2*dvydt3 + dvydt4)/6
        
        return tt, ry, vy
    