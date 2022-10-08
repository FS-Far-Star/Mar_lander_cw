# uncomment the next line if running in a notebook
# %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

# mass, spring constant, initial position and velocity
m = 1
k = 1
x = 0
v = 1

# simulation time, timestep and time
t_max = 100
dt = 0.1                                  # Waves of different frequencies will form if dt>1 for Verlet        
t_array = np.arange(0, t_max, dt)
x_verlet,v_verlet = [x],[v]

# initialise empty lists to record trajectories
x_list = []
v_list = []

# Euler integration
for t in t_array:
    # append current state to trajectories
    x_list.append(x)
    v_list.append(v)

    # calculate new position and velocity
    a = -k * x / m
    x = x + dt * v
    v = v + dt * a

# Borrow the first two terms from Euler
x_verlet.append(x_list[1])
v_verlet.append(v_list[1])

# Verlet Method
for i in range(1,len(x_list)-1):
    a = -k * x_verlet[i] / m
    # x_verlet[i+1] = x_verlet[i] + v_verlet[i] * dt + (a*(dt**2)*0.5)
    # v_verlet[i+1] = v_verlet[i] + a * dt
    x = 2 * x_verlet[i] - x_verlet[i-1] + a*(dt**2)
    v = (x - x_verlet[i-1]) / (2*dt)
    x_verlet.append(x)
    v_verlet.append(v)
print(x_verlet)

# convert trajectory lists into arrays, so they can be sliced (useful for Assignment 2)
x_array = np.array(x_list)
v_array = np.array(v_list)
x_array2 = np.array(x_verlet)
v_array2 = np.array(v_verlet)

# plot the position-time graph
plt.figure(1)
plt.clf()
plt.title("Euler Method")
plt.xlabel('time (s)')
plt.grid()
plt.plot(t_array, x_array, label='x (m)')
plt.plot(t_array, v_array, label='v (m/s)')
plt.legend()

plt.figure(2)
plt.clf()
plt.title("Verlet Method")
plt.xlabel('time (s)')
plt.grid()
plt.plot(t_array, x_array2, label='x (m)')
plt.plot(t_array, v_array2, label='v (m/s)')
plt.legend()

plt.show()
