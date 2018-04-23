#!/apps/anaconda2/bin/python
import matplotlib.pyplot as plt
import numpy as np
import math

"""
Python 2 
Solves 1D heat equation:
	dT/dt = mu(T)*(d^T/dx^2)
	where mu(T) = 1 + 0.001T^Pi

	0 <= x <= 1

Boundary conditions:
	T(t,0) = 1 # Left side
	T(t,1) =  0 # Right


"""


# Stability criteria for Forward-Time Central Space (FTCS) is CFL <= 0.5
cfl = 0.05 

# Space variables
x = 1.0
nx = 51
dx = x/(nx-1)

# Temporal variables
final_time = 0.02
dt = (cfl * dx * dx)/(1.001**math.pi)
nt = int((final_time/dt) - 1)

cfl_new = dt/(dx * dx)  # re-compute CFL number using space and time only

# 1D solution matrix
sol = np.zeros(nx) # time and space grid

print "dx = {0} dt = {1}".format(dx, dt)

# Set initial condition
sol[:] = 0 # T(0,x) = 0      i.e T = 0 @ t = 0, at all x locations

# Set boundary conditions
sol[0] = 1  	# left boundary condition T(t,0) = 1
sol[-1] = 0  	# right boundary condition T(t,1) = 0 Note x: 0 -> 1

for t in range(0, nt-1): # time loop
	for i in range(1, nx-1): # space loop
		sol[i] = sol[i] + (1 + 0.001 * math.pow(sol[i], math.pi)) * cfl_new * (sol[i+1] - 2*sol[i] + sol[i-1])


for i in range(0, nx):
	print sol[i]

l = np.linspace(0, 1, nx)


plt.plot(l, sol[:], linewidth=2, linestyle='--')
plt.xlabel('Length', fontsize=14)
plt.ylabel('Temperature (T)', fontsize=14)
plt.title('1D Heat Equation  @ t = {0} Seconds'.format(final_time), fontsize=18)
plt.grid()
plt.show()

