import numpy as np # maths
import matplotlib.pyplot as plt # plotting

#%matplotlib inline # Uncomment if using Jupyter!


# ========================= Constants & parameters ========================== #

# Constants
g_0 = 9.81				# gravitational acceleration [m/s^2]
rho_0 = 1.225			# air density at sea level [kg/m^3]
H = 8500				# scale height [m]
c_D = 0.54				# drag coefficient for the rocket [-], taken from OpenRocket
A = 0.0083323			# rocket body frontal area [m^2], calculated using the maximum rocket diameter
m_w = 19.803			# “wet” mass of the rocket [kg], data given
m_p = 8.496				# propellant mass [kg], from engine manufacturer
T = 2501.8  			# average rocket engine thrust [N], from engine manufacturer
t_burn = 6.09			# burn time [s], from engine manufacturer

# Simulation parameters
dt = 0.001				# simulation time step [s]
t_f = 150				# simulation time end [s]

# Calculated constants
m_d = m_w - m_p		# “dry” mass of the rocket [kg]
m = m_d + m_p*.5    # average rocket mass [kg]
K = .5*c_D*A/m      # constant in the expression for acceleration due to drag
a_T = T/m           # average acceleration due to thrust


# ================================ Functions ================================ #


def air_density(y):
    """
    Air density [kg/m^3]
    as a function of altitude y [m]
    """
    return rho_0*np.exp(-y/H)


def a_thrust(t):
    """
    Acceleration due to the rocket engine thrust [m/s^2]
    as a function of time t [s]
    """
    return np.where(t < t_burn, a_T, 0)
    # The above line allows our function to be used on entire numpy arrays at
    # once, instead of being restricted to single values. It corresponds to the
    # following if-else statement:
    # if t < t_burn:
    #     a = a_T
    # else:
    #     a = 0
    # return a


def a_drag(y, v):
    """
    Acceleration due to air resistance [m/s^2]
    as a function of altitude y [m] and velocity v [m/s]
    """
    C = K * air_density(y) * np.sqrt(v**2)
    D = - C * v
    return D


def a_gravity(t):
    """
    Acceleration due to gravity [m/s^2]
    as a function of time t [s].
    Note that we assume that the rocket starts at the ground, such that gravity
    is counterbalanced by a normal force at t=0.
    """
    return np.where(t > 0, -g_0, 0)
    # The above line allows our function to be used on entire numpy arrays at
    # once, instead of being restricted to single values. It corresponds to the
    # following if-else statement:
    # if t > 0:
    #     a = -g_0
    # else:
    #     a = 0
    # return a


def a_total(t, y, v):
    """
    Total acceleration [m/s^2]
    as a function of time t [s], altitude y [m], and velocity v [m/s]
    """
    return a_thrust(t) + a_drag(y, v) + a_gravity(t)

def q(y, v):
    q = 0.5 * air_density(y) * v **2
    return q

# ======================== Numerical implementation ========================= #

# Calculate the number of data points in our simulation
datapoints = int(np.ceil(t_f/dt))

# Create data lists
# Except for the time list, all lists are initialized as lists of zeros.
t = np.arange(t_f, step=dt) # runs from 0 to t_f with step length dt
y = np.zeros(datapoints)
v = np.zeros(datapoints)


# We will use while loops to iterate over our data lists. For this, we will
# the auxillary variable i to keep track of which element we're looking at.
# The data points are numbered from 0 to datapoints-1
i = 0
i_max = datapoints - 1


# We iterate until the rocket has crashed, or until we reach the lists' end:
while y[i] >= 0 and i < i_max:
    # Position
    y[i+1] = y[i] + v[i] * dt
    
    # Velocity
    v[i+1] = v[i] + a_total(t[i], y[i], v[i]) * dt
    
    # Advance i with 1
    i += 1

dynamic_pressures = []
for j in range(i):
    y1 = y[j]
    v1 = v[j]
    q1 = q(y1, v1)
    dynamic_pressures.append(q1)

# When we exit the loop above, the rocket has crashed.
# Since we don't need the data after i, we redefine our lists to only include
# the points from 0 to i:
t = t[:i]
y = y[:i]
v = v[:i]


# ============================== Data analysis ============================== #
# exercise a
# Maximal values
apogee = np.max(y)
v_burnout = v[int(t_burn/dt)]
h_burnout = y[int(t_burn/dt)]
t_apogee = 0
t_total = 0
for j in range(i):
    if -0.002 < v[j] < 0.002 and t[j] > 0:
        t_apogee = t[j]
        print(t_apogee)
    if -1 < y[j] < 1 and t[j] > 0:
        t_total = t[j]

max_velocity = max(v) #c

times = np.linspace(0,150,1000)
accs = []
ads = []
for t0 in times:
    y1 = y[int(t0*i/1000)]
    v1 = v[int(t0*i/1000)]
    aa = a_total(t0, y1, v1)
    ad = a_drag(y1, v1)
    accs.append(aa)
    ads.append(ad)

plt.plot(t, v)
plt.plot(times, accs)
plt.plot(times, ads)
plt.show()

#exercise d
max_q = max(dynamic_pressures)

n = 0
for q1 in dynamic_pressures:
    if q1 == max_q:
        v_maxq = v[n]
        h_maxq = y[n]
        t_maxq = t[n]
    n += 1
# =========================== Printing of results =========================== #

#exercise a

print('') # Prints a new line. Not necessary, but used to increase readability
print('Apogee: ', round(apogee), 'm')
print('Burnout velocity: ', round(v_burnout,2), 'm/s')
print('Burnout altitude: ', round(h_burnout), 'm')
print('Time at apogee: ', round(t_apogee,2), 's')
print('Total flight-time: ', round(t_total,2), 's')

print('') # Prints a new line. Not necessary, but used to increase readability

#exercise c
print('Maximum velocity: ', round(max_velocity,2), 'm/s', '(This is equal to the burnout velocity)')

#exercise d
print('Maximum dynamic pressure: ', round(max_q,2), 'Pa')
print('Time, altitude and velocity at maximum dynamic pressure: ', round(t_maxq,2), 's', round(h_maxq), 'm', round(v_maxq,2), 'm/s')

# =========================== Plotting of results =========================== #

# Close all currently open figures, so we avoid mixing up old and new figures.
plt.close('all')

# Trajectory
plt.figure()
plt.plot(t, y)
plt.xlabel("Time [s]")
plt.ylabel("Altitude [m]")
#plt.gca().set_aspect(1) # forces a 1:1 aspect ratio
plt.grid(linestyle='--')
plt.show()