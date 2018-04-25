#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import math


"""
Python 2 

THIS IS A WORK IN PROGRESS

"""


def construct_conservative_variables(PV, nx, ny):
	for i in range(0, nx): # space x loop
		for j in range(0, ny): # space y loop
			CSV[i,j,0] = PV[i,j,0]
			CSV[i,j,1] = PV[i,j,0] * PV[i,j,1]
			CSV[i,j,2] = PV[i,j,0] * PV[i,j,2]
			CSV[i,j,3] = PV[i,j,3]/(gamma-1) + 0.5*PV[i,j,0]*(PV[i,j,1]**2+PV[i,j,2]**2)



def construct_fluxes(PV,CSV,nx,ny):
	for i in range(0, nx): # space x loop
		for j in range(0, ny): # space y loop
			FX[i,j,0] = CSV[i,j,1]
			FX[i,j,1] = ((CSV[i,j,1]**2)/PV[i,j,0]) + PV[i,j,3]
			FX[i,j,2] = CSV[i,j,1] * PV[i,j,2]#
			FX[i,j,3] = (CSV[i,j,3] + PV[i,j,3]) * PV[i,j,1]

			FY[i,j,0] = CSV[i,j,2]
			FY[i,j,1] = FX[i,j,2]
			FY[i,j,2] = ((CSV[i,j,2]**2)/PV[i,j,0]) + PV[i,j,3]#
			FY[i,j,3] = (CSV[i,j,3] + PV[i,j,3]) * PV[i,j,2]#


def re_construct_primitive_variables(CSV,nx,ny):
	for i in range(0, nx): # space x loop
		for j in range(0, ny): # space y loop
			PV[i,j,0] = CSV[i,j,0]
			PV[i,j,1] = CSV[i,j,1]/CSV[i,j,0]
			PV[i,j,2] = CSV[i,j,2]/CSV[i,j,0] #
			PV[i,j,3] = PV[i,j,0] * PV[i,j,1] * PV[i,j,1] 
	
	


# Stability criteria for Forward-Time Central Space (FTCS) is CFL <= 0.25 for uniform x,y grid
cfl = 0.05

# Space variables
x = 0.1
y = 0.1
nx = 210
ny = 210
dx = x/(nx-1)
dy = y/(ny-1)

# Temporal variables
final_time = 0.0005
dt = 2 * cfl * dx
nt = int((final_time/dt) - 1)


PV = np.zeros((nx, ny, 4))
CSV = np.zeros((nx, ny, 4))
FX = np.zeros((nx, ny, 4))
FY = np.zeros((nx, ny, 4))


# Set Initial condition ***************************************************

T = 303 # temperature in kelvin
gamma = 1.4
R = 287

PV[:,:,0] = 1.01
PV[:,:,1] = 0 # u = 0 m/s
PV[:,:,2] = 0
#PV[:,:,3] = PV[:,:,0] * R * T  # density * gas constant * temperature
PV[:,:,3] = PV[:,:,0] * PV[:,:,1] * PV[:,:,1] 

construct_conservative_variables(PV,nx,ny)
construct_fluxes(PV,CSV,nx,ny)

# End of set initial conditions *******************************************


# Set boundary condition **************************************************

# Left
#PV[:,:,0] = 1.01
#PV[0,:,1] = 0 # u = 1 m/s
PV[0,:,2] = 0
#PV[0,:,3] = PV[0,:,0] * PV[0,:,1] * PV[0,:,1] 

# Right
#PV[-1,:,1] = 0 # u = 1 m/s
PV[-1,:,2] = 0
#PV[-1,:,3] = PV[-1,:,0] * PV[-1,:,1] * PV[-1,:,1] 

# Top
PV[:,-1,1] = 1 # u = 1 m/s
#PV[:,-1,2] = 0
PV[:,-1,3] = PV[:,-1,0] * PV[:,-1,1] * PV[:,-1,1] 

# Bottom
PV[:,0,1] = 0 # u = 1 m/s
#PV[:,0,2] = 0
#PV[:,0,3] = PV[:,0,0] * PV[:,0,1] * PV[:,0,1] 

construct_conservative_variables(PV,nx,ny)
construct_fluxes(PV,CSV,nx,ny)


# *************************************************************************


print "Solution Time = {0} dt = {1} dx = {2}  dy = {3}".format(final_time,dt, dx, dy)
print "nt = {0} nx = {1}  ny = {2}  nt*nx*ny = {3}".format(nt, nx,ny, nt*nx*ny)



for t in range(0, nt-1): # time loop
	for i in range(1, nx-1): # space x loop
		for j in range(1, ny-1): # space y loop
			CSV[i,j,:] = CSV[i,j,:] + cfl * ((FX[i-1,j,:] - FX[i+1,j,:]) + (FY[i,j-1,:] - FY[i,j+1,:]))

	re_construct_primitive_variables(CSV,nx,ny)	

	construct_fluxes(PV,CSV,nx,ny)
			

# for plotting
lx = np.linspace(0, x, nx)
ly = np.linspace(0, y, ny)

Y, X = np.meshgrid(ly, lx)


contours = plt.contourf(X, Y, PV[:,:,1],cmap='jet')
plt.colorbar(contours)
plt.clabel(contours, inline=True, fontsize=12, colors='black')
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.title('2D Euler Equation @ t = {} Seconds'.format(final_time), fontsize=18)
plt.show()




