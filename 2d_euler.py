#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import math

"""
Python 2 

*** THIS IS A WORK IN PROGRESS ****

"""

def set_boundary_conditions(t):
    # Primitive variables ******************************************************

	#@ Left
	PV[t,0,:,0] = 1.01 
	#PV[t,0,:,1] = 0 # u = 100 m/s, update all other affected variables and fluxes
	PV[t,0,:,2] = 0
	PV[t,0,:,3] = 100000

	CSV[t,0,:,0] = PV[t,0,:,0] 
	CSV[t,0,:,1] = PV[t,0,:,0] * PV[t,0,:,1]
	CSV[t,0,:,2] = PV[t,0,:,0] * PV[t,0,:,2]
	CSV[t,0,:,3] = PV[t,0,:,3]/(gamma-1) + 0.5*PV[t,0,:,0]*(PV[t,0,:,1]**2+PV[t,0,:,2]**2)


	#@ Bottom
	PV[t,:,0,0] = 1.01
	PV[t,:,0,1] = 0
	#PV[t,:,0,2] = 0
	PV[t,:,0,3] = 100000

	CSV[t,:,0,0] = PV[t,:,0,0] 
	CSV[t,:,0,1] = PV[t,:,0,0] * PV[t,:,0,1]
	CSV[t,:,0,2] = PV[t,:,0,0] * PV[t,:,0,2]
	CSV[t,:,0,3] = PV[t,:,0,3]/(gamma-1) + 0.5*PV[t,:,0,0]*(PV[t,:,0,1]**2+PV[t,:,0,2]**2)

	#@ Top
	PV[t,:,-1,0] = 1.01
	PV[t,:,-1,1] = 1
	#PV[t,:,-1,2] = 0
	PV[t,:,-1,3] = 100000

	CSV[t,:,-1,0] = PV[t,:,-1,0] 
	CSV[t,:,-1,1] = PV[t,:,-1,0] * PV[t,:,-1,1]
	CSV[t,:,-1,2] = PV[t,:,-1,0] * PV[t,:,-1,2]
	CSV[t,:,-1,3] = PV[t,:,-1,3]/(gamma-1) + 0.5*PV[t,:,-1,0]*(PV[t,:,-1,1]**2+PV[t,:,-1,2]**2)

	#@ Right (Exit)
	PV[t,-1,:,0] = 1.01 
	#PV[t,-1,:,1] = 0 # u = 100 m/s, update all other affected variables and fluxes
	PV[t,-1,:,2] = 0
	PV[t,-1,:,3] = 100000

	CSV[t,-1,:,0] = PV[t,-1,:,0] 
	CSV[t,-1,:,1] = PV[t,-1,:,0] * PV[t,-1,:,1]
	CSV[t,-1,:,2] = PV[t,-1,:,0] * PV[t,-1,:,2]
	CSV[t,-1,:,3] = PV[t,-1,:,3]/(gamma-1) + 0.5*PV[t,-1,:,0]*(PV[t,-1,:,1]**2+PV[t,-1,:,2]**2)


def set_boundary_conditions_2():
    # Primitive variables ******************************************************

	#@ Left
	PV[:,0,:,0] = 1.01 
	#PV[:,0,:,1] = 0 # u = 100 m/s, update all other affected variables and fluxes
	PV[:,0,:,2] = 0
	PV[:,0,:,3] = 100000

	CSV[:,0,:,0] = PV[:,0,:,0] 
	CSV[:,0,:,1] = PV[:,0,:,0] * PV[:,0,:,1]
	CSV[:,0,:,2] = PV[:,0,:,0] * PV[:,0,:,2]
	CSV[:,0,:,3] = PV[:,0,:,3]/(gamma-1) + 0.5*PV[:,0,:,0]*(PV[:,0,:,1]**2+PV[:,0,:,2]**2)


	#@ Bottom
	PV[::,0,0] = 1.01
	PV[:,:,0,1] = 0
	#PV[:,:,0,2] = 0
	PV[:,:,0,3] = 100000

	CSV[:,:,0,0] = PV[:,:,0,0] 
	CSV[:,:,0,1] = PV[:,:,0,0] * PV[:,:,0,1]
	CSV[:,:,0,2] = PV[:,:,0,0] * PV[:,:,0,2]
	CSV[:,:,0,3] = PV[:,:,0,3]/(gamma-1) + 0.5*PV[:,:,0,0]*(PV[:,:,0,1]**2+PV[:,:,0,2]**2)

	#@ Top
	PV[:,:,-1,0] = 1.01
	PV[:,:,-1,1] = 1
	#PV[:,:,-1,2] = 0
	PV[:,:,-1,3] = 100000

	CSV[:,:,-1,0] = PV[:,:,-1,0] 
	CSV[:,:,-1,1] = PV[:,:,-1,0] * PV[:,:,-1,1]
	CSV[:,:,-1,2] = PV[:,:,-1,0] * PV[:,:,-1,2]
	CSV[:,:,-1,3] = PV[:,:,-1,3]/(gamma-1) + 0.5*PV[:,:,-1,0]*(PV[:,:,-1,1]**2+PV[:,:,-1,2]**2)

	#@ Right (Exit)
	PV[:,-1,:,0] = 1.01 
	#PV[:,-1,:,1] = 0 # u = 100 m/s, upda:e all o:her affec:ed variables and fluxes
	PV[:,-1,:,2] = 0
	PV[:,-1,:,3] = 100000

	CSV[:,-1,:,0] = PV[:,-1,:,0] 
	CSV[:,-1,:,1] = PV[:,-1,:,0] * PV[:,-1,:,1]
	CSV[:,-1,:,2] = PV[:,-1,:,0] * PV[:,-1,:,2]
	CSV[:,-1,:,3] = PV[:,-1,:,3]/(gamma-1) + 0.5*PV[:,-1,:,0]*(PV[:,-1,:,1]**2+PV[:,-1,:,2]**2)



# Stability criteria for Forward-Time Central Space (FTCS) is CFL <= 0.25 for uniform x,y grid
cfl = 0.1

# Space variables
x = 1.0
y = 1.0
nx = 100
ny = 100
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
PV[0,:,:,1] = 0# u = 100 m/s
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


set_boundary_conditions_2()

#cfl = 0.2

for t in range(0, nt-1): # time loop
	for j in range(1, ny-1): # space y loop
		for i in range(1, nx-1): # space x loop
			#set_boundary_conditions(t)
			CSV[t+1,i,j,:] = CSV[t,i,j,:] + cfl * ((FX[t,i-1,j,:] - FX[t,i+1,j,:]) + (FY[t,i,j-1,:] - FY[t,i,j+1,:]))

			'''
			
			# Update primitive variables
			PV[t+1,i,j,0] = CSV[t+1,i,j,0]
			PV[t+1,j,1] = CSV[t+1,i,j,1]/PV[t+1,i,j,0]
			PV[t+1,i,j,2] = CSV[t+1,i,j,2]/PV[t+1,i,j,0]
			PV[t+1,i,j,3] = PV[t+1,i,j,0] * R * T  # density * gas constant * temperature


			FX[t+1,i,j,0] = CSV[t+1,i,j,1]
			FX[t+1,i,j,1] = CSV[t+1,i,j,1] * PV[t+1,i,j,1] + PV[t+1,i,j,3]
			FX[t+1,i,j,2] = CSV[t+1,i,j,1] * PV[t+1,i,j,2]
			FX[t+1,i,j,3] = (CSV[t+1,i,j,3] + PV[t+1,i,j,3]) * PV[t+1,i,j,1]

			FY[t+1,i,j,0] = CSV[t+1,i,j,2]
			FY[t+1,i,j,1] = CSV[t+1,i,j,1] * PV[t+1,i,j,2]
			FY[t+1,i,j,2] = CSV[t+1,i,j,2] * PV[t+1,i,j,2] + PV[t+1,i,j,3]
			FY[t+1,i,j,3] = (CSV[t+1,i,j,3] + PV[t+1,i,j,3]) * PV[t+1,i,j,2]
			'''
			#FX[t,i,j,0] = CSV[t,i,j,1]

			#FY[t,i,j,0] = CSV[t,i,j,2]
			#PV[t,i,j,1] = CSV[t,i,j,1]/PV[t,i,j,0]
	
	PV[t+1,:,:,0] = CSV[t+1,:,:,0]
	PV[t+1,:,:,1] = CSV[t+1,:,:,1]/PV[t+1,:,:,0]
	PV[t+1,:,:,2] = CSV[t+1,:,:,2]/PV[t+1,:,:,0]
	PV[t+1,:,:,3] = PV[t+1,:,:,0] * R * T  # density * gas constant * temperature


	FX[t+1,:,:,0] = CSV[t+1,:,:,1]
	FX[t+1,:,:,1] = CSV[t+1,:,:,1] * PV[t+1,:,:,1] + PV[t+1,:,:,3]
	FX[t+1,:,:,2] = CSV[t+1,:,:,1] * PV[t+1,:,:,2]
	FX[t+1,:,:,3] = (CSV[t+1,:,:,3] + PV[t+1,:,:,3]) * PV[t+1,:,:,1]

	FY[t+1,:,:,0] = CSV[t+1,:,:,2]
	FY[t+1,:,:,1] = CSV[t+1,:,:,1] * PV[t+1,:,:,2]
	FY[t+1,:,:,2] = CSV[t+1,:,:,2] * PV[t+1,:,:,2] + PV[t+1,:,:,3]
	FY[t+1,:,:,3] = (CSV[t+1,:,:,3] + PV[t+1,:,:,3]) * PV[t+1,:,:,2]
	
print "dt = {0} dx = {1}  dy = {2}".format(dt, dx, dy)



# for plotting
lx = np.linspace(0, x, nx)
ly = np.linspace(0, y, ny)

Y, X = np.meshgrid(ly, lx)


contours = plt.contourf(X, Y, CSV[-1,:,:,1],10,cmap='jet')
plt.colorbar(contours)
plt.clabel(contours, inline=True, fontsize=12, colors='black')
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.title('2D Euler @ t = {} Seconds'.format(tot_time), fontsize=18)
plt.show()




