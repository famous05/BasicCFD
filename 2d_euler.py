#!/apps/anaconda2/bin/python
import matplotlib.pyplot as plt
import numpy as np
import math

"""
Python 2 

THIS IS A WORK IN PROGRESS

"""


# Stability criteria for Forward-Time Central Space (FTCS) is CFL <= 0.25 for uniform x,y grid
cfl = 0.1

# Space variables
x = 10.0
y = 1.0
nx = 51
ny = 51
dx = x/(nx-1)
dy = y/(ny-1)


PV = np.zeros((nt, nx, ny, 4))
CSV = np.zeros((nt, nx, ny, 4))
FX = np.zeros((nt, nx, ny, 4))
FY = np.zeros((nt, nx, ny, 4))


# Set Initial condition ***************************************************

T = 303 # temperature in kelvin
gamma = 1.4
R = 287

PV[0,:,:,0] = 1.01 # density
PV[0,:,:,1] = 0 # u = 0 m/s
PV[0,:,:,2] = 0
PV[0,:,:,3] = PV[0,:,:,0] * R * T  # density * gas constant * temperature

P = PV[0,:,:,3]  # initial pressure
rho  = PV[0,:,:,0]
u = PV[0,:,:,1]
v = PV[0,:,:,2]

CSV[0,:,:,0] = PV[0,:,:,0]
CSV[0,:,:,1] = PV[0,:,:,0] * PV[0,:,:,1]
CSV[0,:,:,2] = PV[0,:,:,0] * PV[0,:,:,2]
CSV[0,:,:,3] = P/(gamma-1) + 0.5*rho*(u*u+v*v)






# Temporal variables
final_time = 0.02
dt = (cfl * dx * dx)/(1.001**math.pi) # assume uniform spacing dx = dy
nt = int((final_time/dt) - 1)

cfl_new = dt/(dx * dx)  # re-compute CFL number using space and time only

# 3D solution matrix with time as 'z'
sol = np.zeros((nt,ny,nx)) # time and space x-y grid

print "dt = {0} dx = {1}  dy = {2}".format(dt, dx, dy)

# Set initial condition
sol[0,:,:] = 0 # T(0,x,y) = 0      i.e T = 0 @ t = 0, at all x and y locations

# Set bounary conditions
sol[:,:,0] = 1 # left side boundary condition T(t,j,0) = 1
sol[:,:,-1] = 0 # right side  boundary condition T(t,j,1) = 0 Note x: 0 -> 1, y: 0-> 1
sol[:,0,:] = 1 # bottom
sol[:,-1,:] = 0 # top


for t in range(0, nt-1): # time loop
	for j in range(1, ny-1): # space y loop
		for i in range(1, nx-1): # space x loop
			sol[t+1,j,i] = sol[t,j,i] + (1 + 0.001 * math.pow(sol[t,j,i], math.pi)) * cfl_new * ((sol[t,j,i+1] - 2*sol[t,j,i] + sol[t,j,i-1]) + (sol[t,j+1,i] - 2*sol[t,j,i] + sol[t,j-1,i]))


# for plotting
lx = np.linspace(0, 10, nx)
ly = np.linspace(0, 1, ny)

Y, X = np.meshgrid(ly, lx)


contours = plt.contourf(Y, X, sol[-1,:,:],11,cmap='jet')
plt.colorbar(contours)
plt.clabel(contours, inline=True, fontsize=12, colors='black')
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.title('2D Heat Equation @ t = {} Seconds'.format(final_time), fontsize=18)
plt.show()




