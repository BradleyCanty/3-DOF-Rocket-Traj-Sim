# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex12.py

Description:
    Make the visualizations more informative by showing the different phases of 
    flight of the rocket. Specifically, for a rocket with n stages the phases of 
    flight are given by...
    
    Stage 1:
        - Vertical flight:                  [t_start_abs,t_start_abs + t_steer_start_rel]
        - Pitch over                        [t_start_abs,t_start_abs + t_steer_start_rel,t_start_abs,t_start_abs + t_steer_start_rel + steer_duration_sec]
        - Gravity turn                      [t_start_abs + t_steer_start_rel + steer_duration_sec, t_bo_abs]
    
    For i in 2 to n...
        Stage i:
            - S(i-1) ECO, stage separation: [t_start_abs,t_start_abs + t_coast_duation_sec]
            - ith stage burn:               [t_start_abs + t_coast_duation_sec,t_bo_abs]
                                             
        if i == n:
            - Sn ECO, coast:                [t_bo_abs, <FINAL TIME IN TIME ARRAY>]
    
    where 
        ECO = Engine Cutoff
    
    Steps to implement this, 
    1) In the Rocket class, create a dictionary having phase names as keys and 
       empty arrays as a values, which will be filled with indexes corresponding 
       to the start and stop times of the phase within the times array

    2) Once the times array is created, iterate through each time in the times 
       array and record the indices of the start and end times of each phase in 
       a separate 'phase_idx' array
       
    3) In the plots, iterate through the dictionary and create a line for each
       key, with the times being sliced according to the indices, stored as the
       value, and label it using the key (since it's the phase name)
    
    Make the same plots found in exercise 11 but showing these flight phases.

Created by: Bradley Canty, 2026/07/19

TO DO: 

"""

import numpy as np
from Rocket import Rocket, RocketStage

v2_s1 = RocketStage(
    name='stage 1',
    m_empty=3000,
    m_propellant=8840,
    m_dot=129.4,
    v_exhaust=2050,
    stage_diameter_m=1.65,
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 3, #sec, steering start time relative to stage's engine start, If first stage, this is the pitch-over time
    steer_duration_sec = 10, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 45, #deg, target flight path angle at stage burnout
    )

v2 = Rocket(
    {
        'name'               : 'V2',
		'date0'              : '2020-01-01', #initial date
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'frame'              : 'J2000',      #standard for earth-centered propagation
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([0,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [v2_s1],
        'tspan'              : 1000,
        'dt'                 : .5,
        'm_payload'          : 1000,
        'payload_separates'  : False,
        'launch_azimuth_deg' : 90,   #Normally changes the sign on eta, right now does nothing...
	})

mm3_s1 = RocketStage(
    name='stage 1',
    m_empty=2060, #kg
    m_propellant=20950, #kg
    m_dot=392, #kg/s
    v_exhaust=2635, #m/s
    stage_diameter_m=1.68, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 3, #sec, steering start time relative to stage's engine start
    steer_duration_sec = 15, #sec, duration to apply thrust offset by angle eta from the velocity vector; set to -1 to steer for the engine burn time
    target_fpa_at_bo_deg = 45, #deg, target flight path angle at stage burnout
    )

mm3_s2 = RocketStage(
    name='stage 2',
    m_empty=680, #kg
    m_propellant=6430, #kg
    m_dot=94, #kg/s
    v_exhaust=2824, #m/s
    stage_diameter_m=1.32, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 0, #sec, steering start time relative to stage's engine start
    steer_duration_sec = -1, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 28, #deg, target flight path angle at stage burnout
    )

mm3_s3 = RocketStage(
    name='stage 3',
    m_empty=350, #kg
    m_propellant=3310, #kg
    m_dot=51, #kg/s
    v_exhaust=2844, #m/s
    stage_diameter_m=1.32, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 0, #sec, steering start time relative to stage's engine start
    steer_duration_sec = -1, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 20, #deg, target flight path angle at stage burnout
    )

mm3 = Rocket(
    {
        'name'               : 'Minuteman 3',
		'date0'              : '2020-01-01', #initial date
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'frame'              : 'J2000',      #standard for earth-centered propagation
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([0,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [mm3_s1,mm3_s2,mm3_s3],
        'tspan'              : 2*60*60,
        'dt'                 : .5,
        'm_payload'          : 1250, #kg, 750 kg mass of PBV and shroud + 500 kg mass of 1 RV
        'payload_separates'  : True,
        'launch_azimuth_deg' : 90,
	})

f9_s1 = RocketStage(
    name='stage 1',
    m_empty=22200, #kg
    m_propellant=411000, #kg, CHECK THIS
    m_dot=2750, #kg/s
    v_exhaust=2800, #m/s
    stage_diameter_m=3.7, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 5, #sec, steering start time relative to stage's engine start
    steer_duration_sec = 0.5, #sec, duration to apply thrust offset by angle eta from the velocity vector; set to -1 to steer for the engine burn time
    target_fpa_at_bo_deg = 45, #deg, target flight path angle at stage burnout
    )

f9_s2 = RocketStage(
    name='stage 2',
    m_empty=4500, #kg
    m_propellant=111500, #kg
    m_dot=287, #kg/s
    v_exhaust=3413, #m/s
    stage_diameter_m=3.7, #m
    coast_duration_sec = 11, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 0, #sec, steering start time relative to stage's engine start
    steer_duration_sec = -1, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 0, #deg, target flight path angle at stage burnout
    )

f9 = Rocket(
    {
        'name'               : 'Falcon 9',
		'date0'              : '2020-01-01', #initial date
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'frame'              : 'J2000',      #standard for earth-centered propagation
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([0,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [f9_s1,f9_s2],
        'tspan'              : 2*60*60,
        'dt'                 : .5,
        'm_payload'          : 17000, #kg
        'payload_separates'  : True,
        'launch_azimuth_deg' : 90,
	})

v2.generate_all_plots(output_img=True,output_gif=True,gif_frame_count=50)
mm3.generate_all_plots(output_img=True,output_gif=True,gif_frame_count=150)
f9.generate_all_plots(output_img=True,output_gif=True,gif_frame_count=150)

#v2.plot_2d_traj(rel_to_launch_point=False,output_img=True,output_gif=False,gif_frame_count=150)
#v2.plot_mass_vs_time(output_img=True,show_flight_phases=True)
#mm3.plot_mass_vs_time(output_img=True,show_flight_phases=True)
#v2.plot_mass_vs_time(output_img=True,show_flight_phases=True)







