import numpy as np
import matplotlib.pyplot as plt

# Mars mass, Gravitational Constant
M = 6.42e23
G = 6.6743e-11
R = 3389.5*10**3

def integrate(coord,v,t_max,dt):
    t_array = np.arange(0, t_max, dt)
    coord_list = []
    v_list = []
    # Euler integration
    for t in t_array:
        # append current state to trajectories
        coord_list.append(coord)
        v_list.append(v)

        # calculate new position and velocity
        magnitude = np.sqrt(np.dot(coord,coord))
        # print('mag',magnitude)
        a = - (G * M) * coord/magnitude**3
        coord = coord + dt * v
        v = v + dt * a
        # print('acc',a)
        # print('v',v)
        # print('coord',coord)

    coord_verlet = coord_list[0:2]
    v_verlet = v_list[0:2]

    # Verlet Method
    for i in range(1,len(coord_list)-1):
        magnitude = np.sqrt(np.dot(coord_verlet[i],coord_verlet[i]))
        a = - (G * M) * coord_verlet[i]/magnitude**3
        coord = 2 * coord_verlet[i] - coord_verlet[i-1] + a*(dt**2)
        v = (coord - coord_verlet[i-1]) / (2*dt)
        coord_verlet.append(coord)
        v_verlet.append(v)

    coord_list = np.array(coord_list)
    v_list = np.array(v_list)
    coord_list2 = np.array(coord_verlet)
    v_list2 = np.array(v_verlet)
    return coord_list,coord_list2

# Circular orbit
m = 1000
coord = np.array([0,20428e3])                   
v = np.array([np.sqrt(G*M/20428e3),0])      #geostationary
t_max ,dt = 100000,1

c1,c2 = integrate(coord,v,t_max,dt)

plt.figure(3)
plt.title("Circular")
plt.plot(0,0,marker=".")        # origin
plt.plot(c1[:,0], c1[:,1])      # Euler
plt.plot(c2[:,0], c2[:,1])      # Verlet
ax = plt.gca() 
ax.set_aspect(1)
draw_circle = plt.Circle((0,0),R,fill=False)
ax.add_artist(draw_circle)
ax.legend(['Origin','Euler','Verlet','Mars surface'],loc='upper right')
# plt.show()

# Descend
m = 1000
coord = np.array([0,20428e3])                   
v = np.array([0,0])
t_max ,dt = 10000,1
t_array = np.arange(0, t_max, dt)
c1,c2 = integrate(coord,v,t_max,dt)

plt.figure(4)
plt.title("Descend")
plt.plot(t_array, c1[:,1])      # Euler
plt.plot(t_array, c2[:,1])      # Verlet
ax = plt.gca() 
# ax.set_aspect(1)
ax.legend(['Euler','Verlet'],loc='upper right')
# plt.show()

# Hyperbolic orbit
m = 1000
coord = np.array([0,20428e3])                      
v = np.array([np.sqrt(2*G*M/20428e3),0])
t_max ,dt = 100000,1

c1,c2 = integrate(coord,v,t_max,dt)

plt.figure(5)
plt.title("Hyperbolic")
plt.plot(0,0,marker=".")        # origin
plt.plot(c1[:,0], c1[:,1])      # Euler
plt.plot(c2[:,0], c2[:,1])      # Verlet
ax = plt.gca() 
ax.set_aspect(1)
draw_circle = plt.Circle((0,0),R,fill=False)
ax.add_artist(draw_circle)
ax.legend(['Origin','Euler','Verlet','Mars surface'],loc='upper right')
# plt.show()

# Eliptical orbit
m = 1000
coord = np.array([0,20428e3])                   
v = np.array([0.6*np.sqrt(G*M/20428e3),0])
t_max ,dt = 100000,1

c1,c2 = integrate(coord,v,t_max,dt)

plt.figure(6)
plt.title("Eliptical")
plt.plot(0,0,marker=".")        # origin
plt.plot(c1[:,0], c1[:,1])      # Euler
plt.plot(c2[:,0], c2[:,1])      # Verlet
ax = plt.gca() 
ax.set_aspect(1)
draw_circle = plt.Circle((0,0),R,fill=False)
ax.add_artist(draw_circle)
ax.legend(['Origin','Euler','Verlet','Mars surface'],loc='upper right')
plt.show()