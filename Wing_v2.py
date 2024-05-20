# -*- coding: utf-8 -*-
"""
Created on Sun May 19 16:06:30 2024

@author: Mathesh JK
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%   AIRFOIL DATA

#given in the tutorial
def compute_CL(alpha):
    CL = 2*np.pi*np.sin(alpha)
    return CL


#%%  BIOT SAVART LAW

#from tutorial
def Biot_Savart(coord_wing, coordinates_wake):
    
    #control point
    XP = coord_wing[0]
    YP = coord_wing[1]
    ZP = coord_wing[2]
    
    #wake point 1
    X1 = coordinates_wake[0][0]
    Y1 = coordinates_wake[0][1]
    Z1 = coordinates_wake[0][2]
    
    #wake point 2
    X2 = coordinates_wake[1][0]
    Y2= coordinates_wake[1][1]
    Z2 = coordinates_wake[1][2]
    
    
    R1 = np.sqrt((XP-X1)**2+(YP-Y1)**2+(ZP-Z1)**2)
    R2 = np.sqrt((XP-X2)**2+(YP-Y2)**2+(ZP-Z2)**2)
    
    R12X = (YP-Y1)*(ZP-Z2) - (ZP-Z1)*(YP-Y2)
    R12Y = -(XP-X1)*(ZP-Z2) + (ZP-Z1)*(XP-X2)
    R12Z = (XP-X1)*(YP-Y2) - (YP-Y1)*(XP-X2)
    
    R12_SQR = R12X**2 + R12Y**2 + R12Z**2
    
    R01 = (X2-X1)*(XP-X1) + (Y2-Y1)*(YP-Y1) + (Z2-Z1)*(ZP-Z2)
    R02 = (X2-X1)*(XP-X2) + (Y2-Y1)*(YP-Y2) + (Z2-Z1)*(ZP-Z2)
    
    CORE = 0.0001 #to avoid singularities at vortex core
    if (R12_SQR < CORE ** 2):
        R12_SQR = CORE ** 2
        # GAMMA = 0;

    if (R1 < CORE):
        R1 = CORE
        # GAMMA = 0;

    if (R2 < CORE):
        R2 = CORE
        # GAMMA = 0;
    
    K = (1/(4*np.pi*R12_SQR))*(R01/R1 - R02/R2)
    
    U = K*R12X
    V = K*R12Y
    W = K*R12Z
    
    return [U,V,W]
    
#%%   CONTROL POINT DISTRIBUTION


b = 10  # wing span
n = 10 # Total number of points in the distribution

# Generate the cosine distribution for z values from 0 to pi
theta = np.linspace(0, np.pi, n)
cosine_values = np.cos(theta)

# Scale the cosine values to the wing span [-b/2, b/2]
y_values = b/2 * cosine_values

# Define x and z coordinates (constant for wing)
x_values = np.zeros(n)
z_values = np.zeros(n)

# Create the 3D coordinates
coordinates_wing = np.zeros([n,3])

coordinates_wing[:,1] = y_values


# Adjust DPI setting for higher resolution
dpi_setting = 300  # Set DPI to desired value

# Plotting the distribution for visualization
fig = plt.figure(dpi=dpi_setting)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_values, y_values, z_values, c=y_values, cmap='viridis')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

#%%  WAKE VORTICES DISTRIBUTION

vortex_wing_pairs = [] #tuple that has points on wing corresponding to one HS vortex
for i in range(int(n/2)):
    pair = (coordinates_wing[i], coordinates_wing[-i-1]) #first and last points form one HSV ans so on
    vortex_wing_pairs.append(pair)



U0 = 10 #inflow vel
c = 1 #chord

rho = 1.225

alpha = 5 #geometric aoa
wake_steps = 5 #no of steps in wake
wake_distance = 10*b #total wake distance to compute
wake_ds = wake_distance/wake_steps #distance between steps
t = wake_ds/U0 #time taken between steps



#initialize matrices
tvortex_pts = np.zeros([wake_steps, n, 3]) 
u_matrix = np.zeros([n,n])
v_matrix = np.zeros([n,n])
w_matrix = np.zeros([n,n])

V_eff_matrix = np.zeros(n)
alpha_eff_matrix = np.zeros(n)
CL_matrix = np.zeros(n).transpose()

gamma = (np.ones(n)).transpose() #initialize gamma
GammaNew_dist = np.zeros(n).transpose()
         
max_itr = 100
tol = 1e-9   
itr = 0


for i in range(n):
    for j in range(n):
        for k in range(wake_steps):
            #first point at TE (wakepoints, n , coordinates)
            if k == 0:
                tvortex_pts[k,:,0] = 0.75*c*np.cos(alpha) #x coordinate 
                tvortex_pts[k,:,1] = coordinates_wing[:,1] #y coordinate, same as wing distribution
                tvortex_pts[k,:,2] = -0.75*c*np.sin(alpha) #z, doesnt change
                
                Velocity = Biot_Savart(coordinates_wing[i], tvortex_pts[k,:,:])
                u = Velocity[0]
                v = Velocity[1]
                w = Velocity[2]
                
                u_matrix[i,j] += u
                v_matrix[i,j] += v
                w_matrix[i,j] += w
                
                V_i = np.sqrt(u**2 + v**2 + w**2) #induced velocity
                
                V_eff = np.sqrt((U0**2) + (V_i**2)) #effective velocity
                
                alpha_i = np.arctan(-w/U0) #induced aoa
                
                alpha_eff = alpha - alpha_i #effective aoa
                
                V_eff_matrix[j] += V_eff #matrix containing a sum of all effective velocities
                #alpha_eff_matrix[j] += alpha_eff
            
            else:
                #for j in range(wake_steps):
                    
                Velocity = Biot_Savart(coordinates_wing[i], tvortex_pts[k,:,:])
                u = Velocity[0]
                v = Velocity[1]
                w = Velocity[2]
                
                
                #Velocity = Biot_Savart(coordinates_wing[i], tvortex_pts[0,:,:])
                
                u_matrix[i,j] += u
                v_matrix[i,j] += v
                w_matrix[i,j] += w
                
                V_i = np.sqrt(u**2 + v**2 + w**2)
                
                V_eff = np.sqrt((U0**2) + (V_i**2))
                
                alpha_i = np.arctan(-w/U0)
                
                alpha_eff = alpha - alpha_i
                
                V_eff_matrix[j] += V_eff
                #alpha_eff_matrix[j] += alpha_eff
                
#%%       TRY TO PUT EVERYTHING UNDER ITERATING LOOP                
#u_matrix has induced vel for unit circulation, hence multiplied with actual gamma                
U_matrix = U0 + (u_matrix @ gamma) 
V_matrix = v_matrix @ gamma
W_matrix = w_matrix @ gamma
"""
Vel_i_matrix = np.sqrt(U_matrix**2 + V_matrix**2 + W_matrix**2)

V_eff_matrix = np.sqrt(U0**2+Vel_i_matrix**2)
"""



U_wake = np.linspace(V_eff_matrix[int(n/2)], U0, wake_steps)  #wake velocity distribution

#t = wake_distance/U0
for i in range(n):
    alpha_eff_matrix[i] = np.arctan(-W_matrix[i]/U0)
    for k in range(wake_steps):
        
        tvortex_pts[k,i,0] = tvortex_pts[k-1,i,0] + t*U_wake[k]
        tvortex_pts[k,i,1] = coordinates_wing[i,1]
        tvortex_pts[k,:,2] = -0.75*c*np.sin(alpha)               
            
while itr < max_itr:            
    print(itr)
    for j in range(n):
        if j >0:
            alpha_eff = alpha_eff_matrix[j]
            V_eff = V_eff_matrix[j]
            
            CL = compute_CL(alpha_eff)
            
            #print(CL)
            
            CL_matrix[j] = CL
            
            Lift = 0.5*rho*(V_eff**2)*CL*c*(coordinates_wing[j,1] - coordinates_wing[j-1,1])
            gammaNew = Lift/(rho*V_eff)
            
            GammaNew_dist[j] = gammaNew
            
    gamma = 0.75*gamma + 0.25*GammaNew_dist
    itr += 1
            
    if np.all(np.abs(gamma - GammaNew_dist) < tol):
        print("solution converged")
        break

plt.plot(y_values,GammaNew_dist)

fig = plt.figure(dpi=dpi_setting)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(tvortex_pts[:,:,0], tvortex_pts[:,:,1], tvortex_pts[:,:,2])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

