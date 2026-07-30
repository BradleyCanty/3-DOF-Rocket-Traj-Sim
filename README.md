# 3-DOF Rocket Trajectory Simulation Exercises [WIP]
## Summary
This is a step-by-step tutorial for creating a 3-DOF rocket trajectory simulation in Python. It begins by comparing a numerical simulation of 1D rocket motion with the Tsiolkovsky rocket equation, and concludes with a desktop application that takes a user-supplied rocket specification (launch location, launch azimuth, launch time, mass parameters per stage, engine parameters per stage, desired burnout flight path angle per stage, payload mass, etc.) and outputs... 
1) plots showing the rocket's trajectory metrics
2) an interactive 3D visualization of the rocket's trajectory over the Earth

## Learning Objectives
Upon completion, the student should be able to...
1) understand rocket motion and how orbit is achieved
2) understand how rocket motion is modeled in a computer
3) understand the reference frames and coordinate systems used in astrodynamics
4) become familiar with modeling and simulation in Python
5) develop an appreciation for numerical simulation

## Motivation
Regarding the last learning objective, numerical simulation is extremely important since it helps us to understand the physical world. An instructive case study is comparing the Tsiolkovsky rocket equation with the numerical integration of the forces acting on a rocket along its trajectory. The Tsiolkovsky rocket equation is a useful, albeit limited, model of rocket motion given by:

deltaV = v_exhaust * ln(m_i/m_f)
where
v_exhaust = the velocity of the exhaust gas
m_i       = the initial mass of the rocket
m_f       = the final mass of the rocket
   
It models the change in velocity of an impulsive, variable-mass system, such as a rocket. A rocket exhausts mass at some velocity over some time span, and the rocket receives an equal but opposite impulse (mass*velocity) in the same time
span. By definition, impulse per unit time is force, in this case thrust. This model is useful because it provides intuition on how rockets move. It's limited since the motion of the system is confined to a single dimension, and the system is absent of outside forces (such as aerodynamic drag or gravity). Consequently, it has very limited real-world use. To model the motion of a real- world system (e.g. one which has multiple degrees of freedom, varies over time,
and interacts with external forces) we can account for all forces and moments on the system along the three principle axes. For simplicity, suppose we only want to model the translation, not the rotation, of the body (3-DOF, not 6-DOF). Our model is then comprised of three equations: force in the x-direction, force in the y-direction, and force in the z-direction. We then divide these forces by mass to get accelerations on the system (i.e., apply Newton's 2nd law: F = m*a) and then numerically integrate once over a small, but finite, time step to get the velocity and once again to get the position. In obtaining the velocity and position for a particular time step, we have "solved" the system for that time step. The
solution of the current time step is used as input for the next time step, i.e., the velocities and positions of the previous time step are used in the acceleration equations in the next time step, which are then integrated at that time step, and so on and so forth until we reach the final time we are interested in. This implies the need for a starting condition: we must know the position and velocity of the system at the starting time so that we can use them in the acceleration equations to solve the system at the first time step. Solving the system at each time step yields acceleration, velocity, and position associated with each time step after the starting time up to and including the ending time. Importantly, the three aforementioned equations that model accelerations along the principle axes are typically coupled and nonlinear, and cannot be solved through analytical means. So, the only way to solve them is through numerical integration. But the point is that they can be solved, and upon doing so yield accurate representations of the real world (*cough* requires verification and validation *cough*) such as the trajectory of rocket. So, in short, the physical world can be modeled faithfully through the careful accounting of forces and a computer. This is, by and large, the motivation for doing these exercises.
 
## Usage and Installation
All exercises are implemented in python and are available in the GitHub repo.
Each exercise has an associated main file (rocket_traj_sim_ex.py):
* the header 'Description' section contains instructions for completing the exercise
* the body contains the solution to the exercise; as a student it is incumbent on you to attempt a solution on your own: if you get stuck, simply refer to the working code
**[TO DO]Additionally, a PowerPoint presentation and suggested reading are supplied for each exercise.**

To install the dependencies, install them using the 'requirements.txt' file:
1) Open your terminal or command prompt
2) Navigate to the directory containing the file:\
   cd path/to/project
5) Activate a new virtual environment for this project (strongly recommended to avoid messing up your global Python installation)
6) Run the installation command:\
   pip install -r requirements.txt

## Outline of Exercises:
1) Simulate the 1-D motion of V2 rocket in absence of outside forces using Euler numerical integration method. Compare results with those of the Tsiolkovsky rocket equation.
2) Modify exercise 1 to apply constant gravity and a thrust cutoff when propellant is expended
3) Apply altitude-dependent gravity, and solve the system of ODEs using midpoint method and 4th order Runge Kutta method
4) Add in the effects of drag using a constant drag coefficient Cd = 0.125
5) Compute air density as a function of altitude using the isothermal barometric formula. Compare this model to tabulated U.S. Standard Atmosphere data (https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html#gsc.tab=0) to verify that it's suitable for use
6) Simulate the trajectory of the rocket using the altitude-varying air density
7) Fit speed of sound versus altitude data with piecewise continuous lines. Additionally, fit drag coefficient vs Mach number data (derived from V2 wind tunnel tests) with a natural cubic spline.
8) Implement the Mach number-varying drag coefficient in the 1D rocket trajectory simulation using the curve fits from the previous exercise.
9) Extend the rocket trajectory simulation from 1D to 2D, and implement a pitch program to initiate a gravity turn. Additionally, make plots of...
   * altitude versus downrange distance
   * flight path angle vs time
   * velocity vs time
   * acceleration vs time
   * dynamic pressure vs time
   * mass vs time
   * trajectory and limb of earth
10) Starting from the previous exercise and assuming an Earth Centered Inertial frame, apply the velocity of the rocket relative to the Earth's rotating atmosphere. This involves reformulating the flight path angle calculation.
11) Organize the code to instantiate a 'rocket' object for each simulation. Organize the plotting functionality such that plots are generated using method calls on the rocket objects. Additionally, extend the code to handle rockets with multiple stages. Finally, make a plot showing the trajectory of a multi-stage rocket.
12) Split the rocket trajectory flight into phases and show the phases in the plots. This will make it easier to discern what's going on:\
    Stage 1:
    - vertical flight
    - pitch over
    - gravity turn
    
   For i in remaining stages:
   - ith stage burn
   - ith engine cutoff
Update all plots to show the trajectory colored by the flight phases, and label the phases using a legend.
13) Refactor the atmospheric model to use the U.S. Standard Atmosphere 1976 specification, then calculate pressure as a function of height, p(h). Update the density, temperature, and speed of sound calculations to use this more accurate atmospheric model. Additionally, update the rocket thrust equation to account for the change in atmospheric pressure, that is, account for the "pressure thrust".
The current rocket thrust equation only uses mass flow rate, m_dot, and exhaust velocity, v_exhaust:\
T = m_dot * v_exhaust\
A more physically-accurate rocket thrust equation introduces a term to account for the change in atmospheric pressure with altitude:\
T = m_dot * v_exhaust + (p_exhaust - p(h)) * A_exit\
where\
p_exhaust = engine exhaust gas pressure at the engine exit plane (i.e. at the end of the nozzle)
p(h)      = air pressure as a function of altitude h
A_exit    = area at the engine exit plane (i.e. at the end of the nozzle)\
Update the RocketStage class to have the p_exhaust and A_exit properties required for computing thrust, then update the thrust terms in the get_state_dot() function in the numerical_tools.py file
14) Implement the following coordinate transformations:
    * LLA to ECI
    * ECI to LLA
    * LLA to ECEF
    * ECEF to LLA
    
    Implement the launch location by starting in the Latitude, Longitude, Altitude (LLA) frame, converting it to the Earth Centered, Earth Fixed (ECEF) frame, and then converting it to the Earth Centered Inertial (ECI) frame. Additionally, change the figure of Earth from sphere to WGS84 ellipsoid, and update the calculations for ground range to use Vincenty's inverse formula. Finally, create the following plots to show the rocket's flight in the ECEF frame:
    * 'Trajectory Metrics (ECEF frame)'
    * 'Trajectory Across Limb of Earth (ECEF frame)'

15) Extend the simulation from 2 dimensions to 3 dimensions, and plot a rocket trajectory on a 3D globe.
16) Implement a launch azimuth by converting from ECI to a topocentric-horizon frame (e.g. SEZ), applying the azimuth rotation to the velocity vector, and then converting it back to ECI. Finally, show the trajectory of a rocket having a non-zero launch azimuth on a 3D globe.
17) Create a GUI for user input and display of outputs, then package the program into an executable using Pyinstaller.
