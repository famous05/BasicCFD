#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import math

"""
Python 2 
Solves 2D heat equation:
	dT/dt = mu(T)*(d^T/dx^2  d^T/dy^2)
	where mu(T) = 1 + 0.001T^Pi

"""


# Stability criteria for Forward-Time Central Space (FTCS) is CFL <= 0.5
cfl = 0.25 

# Space variables
x = 1.0
y = 1.0
nx = 51
ny = 51
dx = x/(nx-1)
dy = y/(ny-1)

# Temporal variables
final_time = 0.02
dt = (cfl * dx * dx)/(1.001**math.pi) # assume uniform spacing dx = dy
nt = int((final_time/dt) - 1)

cfl_new = dt/(dx * dx)  # re-compute CFL number using space and time only

# 3D solution matrix with time as 'z'
sol = np.zeros((nt,ny,nx)) # time and space x-y grid

print "dt = {0} dx = {1}  dy = {2}".format(dt, dx, dy)

# set initial condition
sol[0,:,:] = 0 # T(0,x,y) = 0      i.e T = 0 @ t = 0, at all x and y locations


for t in range(0, nt-1): # time loop
	for j in range(0, ny-1): # space y loop
		for i in range(0, nx-1): # space x loop
			sol[t,j,0] = 1 # left side boundary condition T(t,0,j) = 1
			sol[t,j,-1] = 0 # right side  boundary condition T(t,1,j) = 0 Note x: 0 -> 1, y: 0-> 1
			sol[t,0,i] = 0 # bottom
			sol[t,-1,i] = 0 # top

			sol[t+1,j,i] = sol[t,j,i] + \
			(1 + 0.001 * math.pow(sol[t,j,i], math.pi)) * cfl_new * ((sol[t,j,i+1] - 2*sol[t,j,i] + sol[t,j,i-1]) + (sol[t,j+1,i] - 2*sol[t,j,i] + sol[t,j-1,i]))



# for plotting
lx = np.linspace(0, 1, nx)
ly = np.linspace(0, 1, ny)

Y, X = np.meshgrid(ly, lx)


contours = plt.contourf(Y, X, sol[-1,:,:],10)
plt.colorbar(contours)
plt.clabel(contours, inline=True, fontsize=12)
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.title('2D Heat Equation @ t = {} Seconds'.format(final_time), fontsize=18)
plt.show()




