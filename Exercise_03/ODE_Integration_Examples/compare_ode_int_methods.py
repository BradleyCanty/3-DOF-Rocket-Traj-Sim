# -*- coding: utf-8 -*-
"""
solve_ode_example.py

Description:
    This compares the following integration methods:
    - Tangential (a.k.a. RK1 or Euler)
    - Midpoint (a.k.a. RK2)
    - Runge-Kutta 4 (a.k.a. RK4)
    
    A proper comparison requires solving a differential equation with a known
    solution so that the error of each method can be computed. 
    
    Additionally, the number of calls to the derivative function are set to be 
    the same for each integration method to make it a fair comparison.
    - Tangential    -> 1 derivative function call per iteration
    - Midpoint      -> 2 derivative function calls per iteration
    - Runge-Kutta 4 -> 4 derivative function calls per iteration
    
    The trial ODE to be solved is: dy/dt = 0.3y
    The exact solution is: y = e^(0.3t)
    
Created by:
    Bradley Canty, 2026/06/24
"""

from numerical_tools import propagate_ode
from matplotlib import pyplot as plt
import numpy as np

f = lambda t,y: .3*y
t0 = 0 
y0 = 1 

t_final = 16 #final time, arbitrarily chosen
n = 16 #number of steps, arbitrarily chosen

#Establish number of derivative calls per step for each integration method
ndcps_euler = 1 
ndcps_midpoint = 2
ndcps_rk4 = 4

#Compute the number of steps for each method (required for fair comparison)
dt_exact = 1 #arbitrarily chosen

n_euler = n/ndcps_euler
dt_euler = (t_final - t0)/n_euler

n_midpoint = n/ndcps_midpoint
dt_midpoint = (t_final - t0)/n_midpoint

n_rk4 = n/ndcps_rk4
dt_rk4 = (t_final - t0)/n_rk4

t = np.arange(t0,t_final+dt_exact/2,dt_exact)
y_exact = lambda t:  np.exp(0.3*t)

t_euler, state_vec_euler = propagate_ode(f, [y0], [t0,t_final], dt_euler, 'euler')
t_midpoint, state_vec_midpoint = propagate_ode(f, [y0], [t0,t_final], dt_midpoint, 'midpoint')
t_rk4, state_vec_rk4 = propagate_ode(f, [y0], [t0,t_final], dt_rk4, 'rk4')

#Extract the solutions from the state vectors (trivial since we pass in a single-element state vector)
y_euler = state_vec_euler[:,0]
y_midpoint = state_vec_midpoint[:,0]
y_rk4 = state_vec_rk4[:,0]

plt.plot(t_euler,y_euler,color='red',marker='o',markerfacecolor='none',markeredgecolor='red',label='euler')
plt.plot(t_midpoint,y_midpoint,color='green',marker='o',markerfacecolor='none',markeredgecolor='green',label='midpoint')
plt.plot(t_rk4,y_rk4,color='blue',marker='o',markerfacecolor='none',markeredgecolor='blue',label='4th order Runge-Kutta')
plt.plot(t,y_exact(t),'k',label='exact solution, $y = e^{0.3t}$')
plt.grid(True)
plt.legend()
plt.title(r'Comparison of Numerical Solutions to $\frac{dy}{dt} = 0.3y$')
plt.xlabel('t')
plt.ylabel('y')

#Compute errors and output to console
y_exact_final = y_exact(t_final)
y_euler_final = y_euler[-1]
euler_perc_error = abs((y_euler_final - y_exact_final)/y_exact_final) * 100
y_midpoint_final = y_midpoint[-1]
midpoint_perc_error = abs((y_midpoint_final - y_exact_final)/y_exact_final) * 100
y_rk4_final = y_rk4[-1]
rk4_perc_error = abs((y_rk4_final - y_exact_final)/y_exact_final) * 100

print('ERRORS:')
print(f'Euler:         {round(euler_perc_error,2)}%')
print(f'Midpoint:      {round(midpoint_perc_error,2)}%')
print(f'Runge-Kutta 4: {round(rk4_perc_error,2)}%')

