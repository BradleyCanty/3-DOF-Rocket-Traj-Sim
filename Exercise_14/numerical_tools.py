# -*- coding: utf-8 -*-
"""
numerical_tools.py
"""

import numpy as np
import ode_tools as ot
import geophysical_tools as gt 
import coordinate_conversion_tools as cct
from symbols import deg_symbol, eta_symbol
import warnings
import sys

def get_state_dot(t,dt,state,stages,get_metrics=False,payload_separates=False):
    '''
    Establish system of ODEs to solve, i.e., the time derivative of the
    state vector
    
    Parameters:
    -----------
    t : float array
        array of times at which to evaluate the function
    state:
        state vector given by (x, y, vx, vy)
    
    Returns:
    -------
    float array
        state vector derivative given by (vx,vy,ax,ay)  
    '''
    rx,ry,vx,vy,m = state
    rz = 0 #TEMPROARY DUMMY VALUE
    
    stage = stages[-1]
    for trial_stage in stages:
        if t >= trial_stage.t_start_abs and t < trial_stage.t_start_abs + trial_stage.coast_duration_sec + trial_stage.engine_burn_duration_sec:
            stage = trial_stage
            break
    
    t_start_abs          = stage.t_start_abs
    t_steer_start_rel    = stage.t_steer_start_rel
    steer_duration_sec   = stage.steer_duration_sec
    coast_duration_sec   = stage.coast_duration_sec
    engine_burn_duration_sec   = stage.engine_burn_duration_sec
    required_eta_deg           = stage.required_eta_deg
    
    r = np.array([rx,ry])
    r_3 = np.array([rx,ry,rz])
    r_norm = np.linalg.norm(r)
    v_inertial = np.array([vx,vy])
    v_inertial_3 = np.array([vx,vy,0])
    omega_earth = np.array([0,0,gt.OMEGA_EARTH])
    v_relative = v_inertial_3 - np.cross(omega_earth,r_3) 
    v_relative = v_relative[:2] #slice applied b/c its 2D
    v_relative_norm = np.linalg.norm(v_relative)
    h = r_norm - cct.eci_to_R(rx,ry,rz)
    
    v_steer = get_blended_vel_vec(h,v_relative,v_inertial)
    v_steer_norm = np.linalg.norm(v_steer)
    fpa_rad = get_flight_path_angle_rad(r,v_steer)
    fpa_deg = fpa_rad * 180/np.pi
    
    eta_deg = 0 #deg, angle measured from velocity vector to thrust vector
    eta_rad = eta_deg * np.pi/180
    
    if t < t_start_abs + coast_duration_sec + t_steer_start_rel:
        #Thrust and drag are initially applied along the position vector (not along velocity vector since
        #magnitude of velocity is initially zero, which would result in division by zero in v_steer/v_steer_norm terms)
        thrust_dir_scaling = r/r_norm
        drag_dir_scaling = r/r_norm
    
    if required_eta_deg == None:
        print(required_eta_deg)
        pass
    
    elif t >= t_start_abs + coast_duration_sec + t_steer_start_rel:
        if t < t_start_abs + coast_duration_sec + t_steer_start_rel + steer_duration_sec:
            eta_deg = required_eta_deg
            eta_rad = eta_deg*np.pi/180
        
        #Thrust is applied along the direction of thrust deflection from the velocity vector
        R = np.array([[np.cos(eta_rad),-np.sin(eta_rad)],[np.sin(eta_rad),np.cos(eta_rad)]])
        thrust_dir_scaling = np.matmul(R,v_steer/v_steer_norm)
        
        #Drag is applied along the relative velocity vector, NOT inertial velocity vector
        drag_dir_scaling = v_relative/v_relative_norm
    
    #Compute gravitational acceleration
    a_grav = -gt.G*gt.m_earth/r_norm**2 * r/r_norm
    
    #Compute thrust acceleration
    if t >= t_start_abs and t < coast_duration_sec:
        m_dot = 0
        F_thrust = 0
        m_total = m
    elif (t >= t_start_abs + coast_duration_sec and t < t_start_abs + coast_duration_sec + engine_burn_duration_sec):
        m_dot = stage.m_dot

        F_thrust = stage.v_exhaust * m_dot + (stage.p_exhaust - gt.get_air_pressure(h)) * stage.A_exit
        m_total = m - m_dot*dt
    
        #if staging occurs within the time step before engine burnout, then account for the loss of mass as m_dot
        if t > t_start_abs + coast_duration_sec + engine_burn_duration_sec - dt and t < t_start_abs + coast_duration_sec + engine_burn_duration_sec:
            if stage == stages[-1] and payload_separates == False:
                m_dot = stage.m_dot
                m_total = stage.m_empty + stage.m_other
            else:
                m_dot = stage.m_empty/dt + stage.m_dot #must divide by dt to get the full amount upon integration
                m_total = stage.m_other
    else:
        m_dot = 0
        F_thrust = 0 
        if stage == stages[-1] and payload_separates == False:
            m_total = stage.m_empty + stage.m_other
        else:
            m_total = stage.m_other
            
    a_thrust = F_thrust / m_total * thrust_dir_scaling
    
    #Compute drag acceleration
    Cd = get_drag_coeff(rx,ry,rz,v_relative_norm)
    a_drag = - (1/2 * gt.get_air_density(h) * v_relative_norm**2 * stage.A * Cd) / m_total * drag_dir_scaling
    
    a = a_grav + a_thrust + a_drag
    
    state_dot = np.array([v_inertial[0],v_inertial[1], a[0], a[1], -m_dot])
    
    if get_metrics:
        M = gt.get_mach_number(v_relative_norm,h)
        dynamic_pressure = .5 * gt.get_air_density(h) * v_relative_norm**2
        metrics = np.array([a_grav[0],a_grav[1],a_thrust[0],a_thrust[1],a_drag[0],a_drag[1],Cd,M,fpa_deg,eta_deg,dynamic_pressure])
        return state_dot,metrics
    else:
        return state_dot

def propagate_traj(state0,tspan,dt,stages,launch_azimuth_deg,payload_separates=False,method='rk4',get_metrics=False):
    ode         = get_state_dot
    func        = ot.methods[method]
    times       = np.arange(0,tspan+dt/2,dt) #dt/2 enables including the end point in the array
    steps       = len(times)
    states      = np.zeros((steps, len(state0)))
    states[0]   = state0 #rx0,ry0,vx0,vy0,m0
    
    solve_pitch_program_etas(func,ode,state0,times,dt,stages,launch_azimuth_deg,payload_separates)
    
    print('SOLVING THE STATE VECTOR FOR EACH TIME STEP... ',end='')
    if get_metrics == False:
        for step in range(steps - 1):
            states[step + 1] = func(ode,times[step], dt, states[step],stages,payload_separates)
            
        rx = states[:,0]
        ry = states[:,1]
        vx = states[:,2]
        vy = states[:,3]
        m  = states[:,4]
            
        return rx, ry, vx, vy, m
    
    else:
        metrics = np.zeros((steps,11))
        rx0 = state0[0]
        ry0 = state0[1]
        rz0 = 0 #TEMPORARY DUMMY VALUE
        vx0 = state0[2]
        vy0 = state0[3]
        vz0 = 0 #TEMPORARY DUMMY VALUE
        m0 =  state0[4]
        
        r0 = np.array([rx0,ry0])
        r0_3 = np.array([rx0,ry0,0])
        r0_norm = np.linalg.norm(r0)
        
        v0_inertial3 = np.array([vx0,vy0,0])
        omega_earth = np.array([0,0,gt.OMEGA_EARTH]) 
        v0_relative = v0_inertial3 - np.cross(omega_earth,r0_3)
        v0_relative_norm = np.linalg.norm(v0_relative)
        
        ax_grav0 = -gt.G*gt.m_earth/r0_norm**2 * rx0/r0_norm
        ay_grav0 = -gt.G*gt.m_earth/r0_norm**2 * ry0/r0_norm

        ax_thrust0 = stages[0].v_exhaust * stages[0].m_dot / m0 * rx0/r0_norm
        ay_thrust0 = stages[0].v_exhaust * stages[0].m_dot / m0 * ry0/r0_norm
        ax_drag0 = 0 
        ay_drag0 = 0 

        Cd0 = get_drag_coeff(rx0,ry0,rz0,v0_relative_norm)
        h0 = r0_norm - cct.eci_to_R(rx0,ry0,rz0)
        M0 = gt.get_mach_number(v0_relative_norm,h0)
        fpa_rad0 = get_flight_path_angle_rad(r0,v0_relative_norm)
        fpa_deg0 = fpa_rad0 * 180/np.pi
        angle_from_vel_vec_to_thrust_vec_deg0 = 0
        dynamic_pressure0 = .5 * gt.get_air_density(h0) * v0_relative_norm**2
        
        metrics[0] = [ax_grav0,ay_grav0,ax_thrust0,ay_thrust0,ax_drag0,ay_drag0,Cd0,M0,fpa_deg0,angle_from_vel_vec_to_thrust_vec_deg0,dynamic_pressure0]
        has_impacted = False
        step_at_impact = 0
        for step in range(steps - 1):
            time = times[step]
            state = states[step]
                
            states[step + 1], metrics[step + 1] = func(ode,time,dt,state,stages,get_metrics,payload_separates)
            
            #Check for impact
            rx_ = states[step+1,0]
            ry_ = states[step+1,1]
            rz_ = 0 #TEMPORARY DUMMY VALUE
            r_ = np.sqrt(rx_**2 + ry_**2)
        
            if (r_ < cct.eci_to_R(rx_,ry_,rz_)):
                step_at_impact = step+1
                has_impacted = True
                break
            
        if has_impacted: 
            #Truncate the states and metrics to the relevant times
            rx = states[:,0][:step_at_impact]
            ry = states[:,1][:step_at_impact]
            vx = states[:,2][:step_at_impact]
            vy = states[:,3][:step_at_impact]
            m  = states[:,4][:step_at_impact]
            
            ax_grav   = metrics[:,0][:step_at_impact]
            ay_grav   = metrics[:,1][:step_at_impact]
            ax_thrust = metrics[:,2][:step_at_impact]
            ay_thrust = metrics[:,3][:step_at_impact]
            ax_drag   = metrics[:,4][:step_at_impact]
            ay_drag   = metrics[:,5][:step_at_impact]
            Cd        = metrics[:,6][:step_at_impact]
            M         = metrics[:,7][:step_at_impact]
            fpa_deg   = metrics[:,8][:step_at_impact]
            eta_deg = metrics[:,9][:step_at_impact]
            dynamic_pressure = metrics[:,10][:step_at_impact]
        else:
            rx = states[:,0]
            ry = states[:,1]
            vx = states[:,2]
            vy = states[:,3]
            m  = states[:,4]
            
            ax_grav   = metrics[:,0]
            ay_grav   = metrics[:,1]
            ax_thrust = metrics[:,2]
            ay_thrust = metrics[:,3]
            ax_drag   = metrics[:,4]
            ay_drag   = metrics[:,5]
            Cd        = metrics[:,6]
            M         = metrics[:,7]
            fpa_deg   = metrics[:,8]
            eta_deg = metrics[:,9]
            dynamic_pressure = metrics[:,10]
        
        print('DONE\n')
        return [rx, ry, vx, vy, m], [ax_grav,ay_grav,ax_thrust,ay_thrust,ax_drag,ay_drag,Cd,M,fpa_deg,eta_deg,dynamic_pressure]


def solve_pitch_program_etas(func,ode,state0,times,dt,stages,launch_azimuth_deg,payload_separates,max_iter=100):
    '''
    Get the thrust steering angle, eta, used in the stage's pitch program
    This goes through each stage and checks whether it has an eta set. If not, the eta is solved for using
    the secant method and the desired burnout flight path angle (found in stage.target_fpa_at_bo_deg)
    '''
    print(f"SOLVING FOR EACH STAGE\'S STEERING ANGLE {eta_symbol}...")
    get_metrics = True
        
    for stage in stages:
        if stage.required_eta_deg == None:
            target_fpa_at_bo_deg = stage.target_fpa_at_bo_deg
            t_pitch_prog_end = stage.t_start_abs + stage.coast_duration_sec + stage.engine_burn_duration_sec
            
            #set the initial eta guess
            prev_eta_guess = 0 #deg
            
            #Set the sign of the second eta guess based on the launch azimuth (this works for 2D simulation only!!!)
            if launch_azimuth_deg == 90:
               curr_eta_guess = 1 #deg
                
            elif launch_azimuth_deg == -90:
                curr_eta_guess = -1
            
            eta_guesses = [prev_eta_guess,curr_eta_guess]
            trial_fpas = [None,None]
            
            iter_ = 1
            while iter_ < max_iter:
                for j,eta_guess in enumerate(eta_guesses):
                    stage.required_eta_deg = eta_guess
                    for i,tt in enumerate(times):
                        if tt >= t_pitch_prog_end:
                            pitch_prog_end_idx = i
                            break
                    
                    times_slice = times[:pitch_prog_end_idx]
                    steps       = len(times_slice)
                    states      = np.zeros((steps, len(state0)))
                    states[0]   = state0 #rx0,ry0,vx0,vy0,m0
                    
                    for step in range(steps - 1):
                        states[step + 1],metrics = func(ode,
                                                        times_slice[step], 
                                                        dt, 
                                                        states[step],
                                                        stages,
                                                        get_metrics,
                                                        payload_separates)
                        

                    trial_fpas[j] = metrics[8]
                
                prev_eta_guess = eta_guesses[0]
                curr_eta_guess = eta_guesses[1]
                prev_trial_fpa = trial_fpas[0]
                curr_trial_fpa = trial_fpas[1]
                
                next_eta_guess = curr_eta_guess + (target_fpa_at_bo_deg - curr_trial_fpa) * (curr_eta_guess - prev_eta_guess)/(curr_trial_fpa-prev_trial_fpa)
                
                stage.required_eta_deg = next_eta_guess
                for i,tt in enumerate(times):
                    if tt >= t_pitch_prog_end:
                        pitch_prog_end_idx = i
                        break
                
                times_slice = times[:pitch_prog_end_idx]
                steps       = len(times_slice)
                states      = np.zeros((steps, len(state0)))
                states[0]   = state0 #rx0,ry0,vx0,vy0,m0
                
                for step in range(steps - 1):
                    states[step + 1],metrics = func(ode,
                                                    times_slice[step], 
                                                    dt, 
                                                    states[step],
                                                    stages,
                                                    get_metrics,
                                                    payload_separates)
                trial_fpa = metrics[8]

                if abs(trial_fpa - target_fpa_at_bo_deg) < .1: #is within 0.1 degrees
                    stage.fpa_at_bo_deg = trial_fpa
                
                    print(f'Solving for {stage.name} steering angle {eta_symbol}, iter = {iter_}')
                    print(f'{eta_symbol} found to be {stage.required_eta_deg}{deg_symbol} in {iter_} iterations')
                    print(f'target_fpa_deg = {target_fpa_at_bo_deg}{deg_symbol}, trial_fpa_deg = {stage.fpa_at_bo_deg}{deg_symbol}\n')
                    break
                
                print(f'Solving for {stage.name} steering angle {eta_symbol}, iter = {iter_}')
                print(f'Previous guessed {eta_symbol} is {prev_eta_guess}{deg_symbol}')
                print(f'Current guessed {eta_symbol} is {curr_eta_guess}{deg_symbol}')
                print(f'Next guessed {eta_symbol} is {next_eta_guess}{deg_symbol}')
                print(f'target_fpa_deg = {target_fpa_at_bo_deg}{deg_symbol}, trial_fpa_deg = {trial_fpa}{deg_symbol}\n')
                
                if abs(next_eta_guess) > 90:
                    print(f'Absolute value of next guessed {eta_symbol} is greater than 90{deg_symbol} but FPA not attained.')
                    print('Failed to converge: try increasing stage steering duration or engine burn time (via increased propellant mass), or change the desired FPA at burnout')
                    sys.exit(1) #end program with exit code 1
                
                prev_eta_guess = curr_eta_guess
                curr_eta_guess = next_eta_guess
                eta_guesses = [prev_eta_guess,curr_eta_guess]
                iter_+=1
                
                if iter_ > max_iter:
                    warn_msg = f'ERROR: Eta solver failed to converge in {max_iter} iterations, exiting program...'
                    warnings.warn(warn_msg)
                    sys.exit(1) #end program with exit code 1

def get_drag_coeff(rx,ry,rz,v_norm):
    r_norm = np.linalg.norm([rx,ry,rz])
    h = r_norm - cct.eci_to_R(rx,ry,rz)
    M = gt.get_mach_number(v_norm, h)
    
    #Data for Cd vs Mach profile of V1 rocket
    M_data = [0.0,0.5,0.7,1,1.15,1.5,1.75,2.0,2.5,3.5,5]
    Cd_data = [.25,.18,.17,.29,.38,.27,.21,.18,.15,.12,.1]
    if M < M_data[-1]:
        return get_cubic_spline_point(M_data,Cd_data,M)
    else: 
        return Cd_data[-1]

def get_blended_vel_vec(h, v_rel, v_inertial, h_min=40000.0, h_max=80000.0):
    """
    Blends the vector smoothly between relative and inertial velocity vectors based on altitude.
    
    Parameters:
        h (float): Current altitude in meters.
        v_rel (np.ndarray): 3D relative velocity vector.
        v_inertial (np.ndarray): 3D inertial velocity vector.
        h_min (float): Altitude where blending begins (meters). Default 40km.
        h_max (float): Altitude where blending ends (meters). Default 80km.
    """
    # 1. Normalize and clamp altitude into a [0, 1] space
    x = (h - h_min) / (h_max - h_min)
    x = max(0.0, min(1.0, x))
    
    # 2. Compute smoothstep weight (Cubic Hermite Spline)
    w = 3 * (x**2) - 2 * (x**3)
    
    # 3. Perform the vector blend
    v_steer = (1.0 - w) * v_rel + w * v_inertial
    
    return v_steer

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

def get_flight_path_angle_rad(r,v):
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    
    #if velocity is zero, then assume FPA is 90 degrees (i.e. on launch pad pointing toward zenith)
    if v_norm == 0:
        return np.pi/2
    
    #Clamp to avoid floating point precision errors outside [-1,1]
    val = np.clip(np.dot(r,v)/(r_norm * v_norm))

    return np.pi/2 - np.acos(val)

def get_radial_vel(r,v):
    array_len = r.shape[0]
    v_r = np.empty(array_len)
    for i in range(array_len):
        r_norm = np.linalg.norm(r[i,:])
        v_r[i] = np.dot(r[i,:],v[i,:])/r_norm
    return v_r

def get_tangential_vel(r,v):
    v_r = get_radial_vel(r,v)
    array_len = r.shape[0]
    v_t = np.empty(array_len)
    for i in range(array_len):
        v_norm = np.linalg.norm(v[i,:])
        if (abs(v_norm**2 - v_r[i]**2) < 10e-5):
            v_t[i] = 0
        else:
            v_t[i] = np.sqrt(v_norm**2 - v_r[i]**2)
    return v_t

def get_radial_accel(r,a):
    array_len = r.shape[0]
    a_r = np.empty(array_len)
    for i in range(array_len):
        r_norm = np.linalg.norm(r[i,:])
        a_r[i] = np.dot(r[i,:],a[i,:])/r_norm
    return a_r

def get_tangential_accel(r,a):
    a_r = get_radial_accel(r,a)
    array_len = r.shape[0]
    a_t = np.empty(array_len)
    for i in range(array_len):
        a_norm = np.linalg.norm(a[i,:])
        if (abs(a_norm**2 - a_r[i]**2) < 10e-5):
            a_t[i] = 0
        else:
            a_t[i] = np.sqrt(a_norm**2 - a_r[i]**2)
    return a_t









