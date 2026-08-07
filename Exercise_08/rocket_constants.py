# -*- coding: utf-8 -*-
"""
rocket_constants.py
"""
import numpy as np

#V1 Rocket
m_empty = 4000 #kg, empty mass (includes 1000kg warhead)
m_propellant = 8840 #kg
m_dot = 129.4 #kg/s, mass flow rate of fuel
v_exhaust = 2050 #m/s, exhaust velocity
D = 1.65 #m, diameter of V2 rocket
A = np.pi*(D/2)**2
#F_thrust = v_exhaust * m_dot
t_burnout = m_propellant/m_dot #seconds
stage_count = 1

pitch_programs = {'t_pitch_prog_start'   :3,         #sec, time when pitch program starts
                  't_pitch_prog_end'     :t_burnout, #sec, time when pitch program ends
                  'steering_duration_sec':10,        #sec, duration to apply thrust offset by angle eta from the velocity vector
                  'target_fpa_at_bo_deg' :30,        #deg, target flight path angle
                  'required_eta_deg'     :None,      #deg, the angle measured from the velocity vector to the thrust vector required to reach the target FPA; to be found
                  'azimuth_deg'          :90}        #deg, measured clockwise from true north
