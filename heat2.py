#importing all necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from math import exp
# define initial conditions
dx=0.2
dt=0.01
x0 = 20
xn = 20
alpha = 1
xmin=0
tmin=0
xmax = 40
tmax = 5
# defining initial conditions for each case
def initial_middle(dx, temp):
    # this one is for having the exact middle become the source of contant temp
    arrx = np.arange(xmin,xmax,dx)
    arry = np.arange(xmin,xmax,dx) # let's have the same dimensions/grid spacing in the "y" direction
    initial_state = np.zeros((len(arrx), len(arry)), dtype=np.float)
    mid_coord = int(len(arrx)/2) 
    initial_state[mid_coord,mid_coord] = temp # set the condition
    return initial_state

def initial_walls(dx,temp):
    # this one is for having the walls become the source of contant temp
    arrx = np.arange(xmin,xmax,dx)
    arry = np.arange(xmin,xmax,dx) # let's have the same dimensions/grid spacing in the "y" direction
    initial_state = np.zeros((len(arrx), len(arry)), dtype=np.float)
    for j in range(len(arrx)):
        for k in range(len(arry)):
            if j==0 or k==0 or j==(len(arrx)-2) or k==(len(arry)-2):
                initial_state[j,k] = temp
    return initial_state

def Heat_2D_next_u(prev_u, stencil, alpha):
    # A generic finite difference method that we use for both boundary conditions,
    # wherein "stencil" is the initial conditions of each boundary condition problem, encoding both the 
    # shape and the values at each place
    u_at_time = np.zeros(prev_u.shape,dtype=np.float) # create a 2D array
    maxx, maxy = prev_u.shape  # set the limits on the iterators
    # we are going to use "j,k" notation to mean the jth row, kth column
    for j in range(maxx-1): # loop through each row
        for k in range(maxy-1): # loop along the columns, down each particular row
            if stencil[j,k] == 20:
                u_j_k_n1 = stencil[j,k] # if the stencil (initial condition) has nonzero place, just use it
                # this works to hold a temp at 20, say, but not at zero. Need to make new mechanism
                # I suppose I could just create a second stencil holding 1s at places with zero...
            else:
                u_j1_k_n = prev_u[j+1,k] # this is the "j+1"th  point of the previous function
                u_j_k1_n = prev_u[j,k+1]
                u_j_1_k_n = prev_u[j-1,k]
                u_j_k_1_n = prev_u[j,k-1]
                u_j_k_n = prev_u[j,k] # previous time step, same location
                # now for the next time step
                u_j_k_n1 = u_j_k_n + alpha*(dt/((dx)**2))*(u_j1_k_n + u_j_k1_n + u_j_1_k_n + u_j_k_1_n - 4*u_j_k_n)
            u_at_time[j,k] = u_j_k_n1
    # now we've looped through each spatial coordinate
    return u_at_time # return the full function after looping through entire spatial interval

# create our "hot center" initial condition
stencil = initial_middle(dx, 20)
u = []  # there will be "n" ndarrays of shape (j,k)
xarrr = np.arange(xmin,xmax,dx) # defines the discretized spatial interval
t = np.arange(tmin,tmax,dt) # defines the discretized temportal interval
for time in t:
    print(time)
    if time == 0.0:
        u.append(stencil) # append the full function to list of functions
    else:
        prev_u = u[-1] # get the previous function
        res = Heat_2D_next_u(prev_u, stencil, alpha)
        u.append(res) # append the full function to list of functions
indice = len(t)-1
plt.imshow(u[indice], cmap='hot', interpolation='nearest', extent=[0,xmax, 0, xmax])
plt.title('2D Heat Diffusion, dx=0.2, dt=0.01')
plt.xlabel('Spatial X-axis')
plt.ylabel('Spatial Y-axis')
plt.show()

# create our "hot walls, cold center" initial condition
stencil2 = initial_walls(dx, 20)
u = []  # remake functions list
xarrr = np.arange(xmin,xmax,dx) # defines the discretized spatial interval
t = np.arange(tmin,tmax,dt) # defines the discretized temportal interval
for time in t:
    print(time)
    if time == 0.0:
        u.append(stencil2) # append the full function to list of functions
    else:
        prev_u = u[-1] # get the previous function
        res = Heat_2D_next_u(prev_u, stencil2, alpha)
        u.append(res) # append the full function to list of functions
indice = len(t)-1
maxt = indice
plt.imshow(u[maxt], cmap='hot', interpolation='nearest', extent=[0,xmax, 0, xmax])
plt.title('2D Heat Diffusion, dx=0.2, dt=0.01')
plt.xlabel('Spatial X-axis')
plt.ylabel('Spatial Y-axis')
plt.show()
