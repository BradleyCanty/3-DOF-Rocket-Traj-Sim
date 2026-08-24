# -*- coding: utf-8 -*-
"""
solve_two_body_orbit.py

This computes a 2-body orbit using numerical integration.
The specific energy of the orbit at the final time is computed and compared with 
the ideal value. The percent error of the specific energy is a good measure of 
the accuracy of the numerical solution.
"""

import numerical_tools as nt
import numpy as np
import matplotlib.pyplot as plt
import os

arrow_symbol = r'$\rightarrow$'

img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)

#EARTH MODEL (WGS84)
G = 6.674E-11 #m^3/(kg*s^2), universal gravitational constant
m_earth = 5.972E24 #kg, mass of earth
R_eq = 6378137.0 # km, Earth equatorial radius
f = 1 / 298.257223563 # flattening
R_p = R_eq*(1-f) # km, Earth polar radius
R_mean = (2*R_eq + R_p)/3 #m, mean radius of Earth


class two_body_orbit():
    
    def __init__(self,name,a,ecc,dt,method):
        self.name = name
        self.a = a
        self.ecc = ecc
        self.dt = dt
        self.method = method
        
        #Compute position vector at perigee
        rx0 = a*(1-ecc) #location at perigee
        ry0 = 0

        #Compute velocity vector at perigee
        vx0 = 0
        vy0 = np.sqrt(G*m_earth/a * (1+ecc)/(1-ecc)) #velocity at perigee, which is normal to the position vector

        #Compute orbital period
        self.p = np.sqrt(4*np.pi**2 * a**3/(G*m_earth))

        self.state0 = [rx0,ry0,vx0,vy0]
        t_final = self.p
        self.tspan = [0,t_final]

    
        self.times,states = nt.propagate_ode(self.state_dot,self.state0,self.tspan,self.dt,method=self.method)
        self.rx = states[:,0]
        self.ry = states[:,1]
        self.vx = states[:,2]
        self.vy = states[:,3]

        #Compare the numerical and theoretical specific energies at end state
        r_norm = np.sqrt(self.rx**2+self.ry**2)
        v_norm = np.sqrt(self.vx**2+self.vy**2)
        self.specific_energy = self.get_specific_energy(r_norm[-1],v_norm[-1])
        self.ideal_specific_energy = self.get_ideal_specific_energy(self.a)
        self.specific_energy_perc_error = abs((self.specific_energy - self.ideal_specific_energy)/self.ideal_specific_energy) * 100    
        
        self.method_name_dict = {'euler'    : 'Euler',
                                 'midpoint' : 'Midpoint',
                                 'rk4'      : 'Runge Kutta 4'}
        
        #Compute altitudes of apogee and perigee
        self.alt_perigee_km = (min(r_norm) - R_mean)/1000
        self.alt_apogee_km = (max(r_norm) - R_mean)/1000
        
    def state_dot(self,t,state):
        r = state[:2]
        accel = -G*m_earth/np.linalg.norm(r)**3 * r
        return np.array([state[2],state[3],accel[0],accel[1]])

    #Compute specific energy
    def get_specific_energy(self,r_norm,v_norm):
        return (v_norm**2)/2 - G * m_earth / r_norm

    def get_ideal_specific_energy(self,a):
        return -G * m_earth / (2*a)
    
    def print_specific_energy_info(self):
        print(f'{self.name} orbit specific energy info:')
        print(f'Numerical solution method: {self.method_name_dict[self.method]} with dt = {self.dt}[s]')
        print(f'Specific Energy at t={round(self.tspan[-1],2)}[s]:')
        print(f'\tIdeal:    {round(self.ideal_specific_energy,2)}[J/kg]')
        print(f'\tComputed: {round(self.specific_energy,2)}[J/kg]')
        print(f'\tPercent error: {round(self.specific_energy_perc_error,6)}%\n')
    
    def plot(self):
        limb_of_earth = plt.Circle((0, 0), R_mean, color='b',zorder=0)
        fig, ax1 = plt.subplots()
        ax1.add_patch(limb_of_earth)
        ax1.plot(self.rx,self.ry,'r')
        ax1.grid(True)
        ax1.axis('scaled')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')

        #Plot times at appropriate intervals along the the trajectory
        time_unit=None
        time_point_interval=None
        time_unit_multiple=None
        if self.times[-1] > 3600 * 4: #if time duration is greater than 4 hours, plot time points in intervals of 1 hour
            time_point_interval = 3600
            time_unit = 'hr'
            time_unit_multiple = 1
        else: #if time duration is less than 4 hours, plot time points in intervals of 10 minutes
            time_point_interval = 600
            time_unit = 'min'
            time_unit_multiple = 10

        for i in range(self.tspan[0],len(self.times),int(time_point_interval/self.dt)):
            ax1.scatter(self.rx[i],self.ry[i],color='black',marker='.',s=100,zorder=2)
            #Orient the text such that it doesn't overlap the orbit
            if self.rx[i] > (max(self.rx) + min(self.rx))/2:
                xtextloc = 30
            elif self.rx[i] == 0:
                xtextloc = 0
            else:
                xtextloc = -30
            if self.ry[i] > (max(self.ry) + min(self.ry))/2:
                ytextloc = 7
            elif self.ry[i] == 0:
                ytextloc = -3
            else:
                ytextloc = -13
                    
            ax1.annotate(
                f't={int(self.times[i]/time_point_interval*time_unit_multiple)}[{time_unit}]', #the text label
                (self.rx[i], self.ry[i]),  # Point coordinate to label
                textcoords="offset points",
                xytext=(xtextloc, ytextloc),
                ha='center'
            )
        p_plot_title=round(self.p/(time_point_interval)*time_unit_multiple,2)
        ax1.set_title(f'{self.name}',fontsize=20,weight='bold')
        plt.xlim(-.25*(max(self.rx) - min(self.rx))+min(self.rx), .25*(max(self.rx) - min(self.rx))+max(self.rx))
        plt.ylim(-.25*(max(self.ry) - min(self.ry))+min(self.ry), .25*(max(self.ry) - min(self.ry))+max(self.ry))
        
        text_str = f"Orbit Parameters:\na={self.a/1000}[km]\ne={self.ecc}\nP={p_plot_title}[{time_unit}]\n$h_{{perigee}}$={round(self.alt_perigee_km,1)}[km]\n$h_{{apogee}}$={round(self.alt_apogee_km,1)}[km]\n\nNumerical Solution Details:\nmethod={self.method_name_dict[self.method]}\ndt={self.dt}[s]\n$E_{{sp}}$ perc. err.={round(self.specific_energy_perc_error,6)}%"
        props = dict(boxstyle='square', facecolor='lightyellow', alpha=0.5)
        
        # Using ax.text with transAxes. X=0.05 (inside left), Y=1.05 (just above the top edge)
        ax1.text(1.12, 0, text_str, transform=ax1.transAxes, fontsize=8, bbox=props, va='bottom')
        
        plt.tight_layout()
        plt.savefig(img_dir + f'/{self.name.replace(' ','_')}_orbit.png', dpi = 600)
        plt.show()
        print(f'Plot image of {self.name} found at:\n{img_dir}\n')

if __name__ == '__main__':
    c1_orbit = two_body_orbit('Circular Orbit RK4',a=R_eq*1.5,ecc=0,dt=100,method='rk4')
    c1_orbit.plot()
    c1_orbit.print_specific_energy_info()
    
    c2_orbit = two_body_orbit('Circular Orbit Euler',a=R_eq*1.5,ecc=0,dt=100,method='euler')
    c2_orbit.plot()
    c2_orbit.print_specific_energy_info()
    
    c2_orbit = two_body_orbit('Circular Orbit Midpoint',a=R_eq*1.5,ecc=0,dt=100,method='midpoint')
    c2_orbit.plot()
    c2_orbit.print_specific_energy_info()
    
    m1_orbit = two_body_orbit('Molniya Orbit RK4',a=26600000,ecc=0.74,dt=100,method='rk4')
    m1_orbit.plot()
    m1_orbit.print_specific_energy_info()
    
    m2_orbit = two_body_orbit('Molniya Orbit Euler',a=26600000,ecc=0.74,dt=100,method='euler')
    m2_orbit.plot()
    m2_orbit.print_specific_energy_info()
    
    m2_orbit = two_body_orbit('Molniya Orbit Midpoint',a=26600000,ecc=0.74,dt=100,method='midpoint')
    m2_orbit.plot()
    m2_orbit.print_specific_energy_info()
    



