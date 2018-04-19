#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import math

"""
Python 2 

THIS IS A WORK IN PROGRESS

"""


# Stability criteria for Forward-Time Central Space (FTCS) is CFL <= 0.25 for uniform x,y grid
cfl = 0.005

# Space variables
x = 10.0
y = 1.0
nx = 1000
ny = 10
dx = x/(nx-1)
dy = y/(ny-1)

tot_time = 0.02
dt = cfl * 2 * dx
nt = int((tot_time/dt) - 1)

print "nt = " +str(nt)

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

FX[0,:,:,0] = CSV[0,:,:,1]
FX[0,:,:,1] = CSV[0,:,:,1] * PV[0,:,:,1] + PV[0,:,:,3]
FX[0,:,:,2] = CSV[0,:,:,1] * PV[0,:,:,2]
FX[0,:,:,3] = (CSV[0,:,:,3] + PV[0,:,:,3]) * PV[0,:,:,1]

FY[0,:,:,0] = CSV[0,:,:,2]
FY[0,:,:,1] = CSV[0,:,:,1] * PV[0,:,:,2]
FY[0,:,:,2] = CSV[0,:,:,2] * PV[0,:,:,2] + PV[0,:,:,3]
FY[0,:,:,3] = (CSV[0,:,:,3] + PV[0,:,:,3]) * PV[0,:,:,2]
# End of set initial condition

# Set boundary condition
#PV[0,0,:,1] = 100 # u = 100 m/s, update all other affected variables and fluxes
#CSV[0,0,:,1] = PV[0,0,:,0] * PV[0,0,:,1]
#CSV[0,0,:,3] = PV[0,0,:,3]/(gamma-1) + 0.5*PV[0,0,:,0]*(PV[0,0,:,1]**2+PV[0,0,:,2]**2)



#cfl = 0.2

for t in range(0, nt-1): # time loop
	for j in range(1, ny-1): # space y loop
		for i in range(1, nx-1): # space x loop
			CSV[t+1,i,j,:] = CSV[t,i,j,:] + cfl * ((FX[t,i-1,j,:] - FX[t,i+1,j,:]) + (FY[t,i,j-1,:] - FY[t,i,j+1,:]))

			#PV[t,i,j,1] = CSV[t,i,j,0]/PV[t,i,j,0]


print "dt = {0} dx = {1}  dy = {2}".format(dt, dx, dy)



# for plotting
lx = np.linspace(0, x, nx)
ly = np.linspace(0, y, ny)

Y, X = np.meshgrid(ly, lx)


contours = plt.contourf(Y, X, PV[-1,:,:,1],11,cmap='jet')
plt.colorbar(contours)
plt.clabel(contours, inline=True, fontsize=12, colors='black')
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.title('2D Euler @ t = {} Seconds'.format(tot_time), fontsize=18)
plt.show()




