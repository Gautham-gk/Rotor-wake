# -*- coding: utf-8 -*-
"""
Created on Thu May 23 10:59:00 2024

@author: Mathesh JK
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



#%%     AIRFOIL DATA

#given in the tutorial
def compute_CL(alpha):
    CL = 2*np.pi*np.sin(alpha)
    return CL


#%%     BIOT SAVART LAW


#from tutorial
def Biot_Savart(coordinates_wing, coordinates_wake_1,coordinates_wake_2):

    #control point
    XP = coordinates_wing[0]
    YP = coordinates_wing[1]
    ZP = coordinates_wing[2]
    #print("Control points are:",XP,YP,ZP)
    
    #wake point 1
    X1 = coordinates_wake_1[0]
    Y1 = coordinates_wake_1[1]
    Z1 = coordinates_wake_1[2]
    #print("Wake point 1 are:",X1,Y1,Z1)
    
    #wake point 2
    X2 = coordinates_wake_2[0]
    Y2= coordinates_wake_2[1]
    Z2 = coordinates_wake_2[2]
    #print("Wake point 2 are:",X2,Y2,Z2)
    
    
    #R1 vector connects wake point 1 and control point
    R1 = np.sqrt((XP-X1)**2+(YP-Y1)**2+(ZP-Z1)**2)

    #R2 vector connects wake point 2 and control point
    R2 = np.sqrt((XP-X2)**2+(YP-Y2)**2+(ZP-Z2)**2)


    #x,y, and z coordinates of R1 X R2
    R12X = (YP-Y1)*(ZP-Z2) - (ZP-Z1)*(YP-Y2)
    R12Y = -(XP-X1)*(ZP-Z2) + (ZP-Z1)*(XP-X2)
    R12Z = (XP-X1)*(YP-Y2) - (YP-Y1)*(XP-X2)

    #Magnitude of R1 X R2
    R12_SQR = R12X**2 + R12Y**2 + R12Z**2


    #Dot poduct of R12(connecting wake points) and R1
    R01 = (X2-X1)*(XP-X1) + (Y2-Y1)*(YP-Y1) + (Z2-Z1)*(ZP-Z1)


    #Dot poduct of R12(connecting wake points) and R2
    R02 = (X2-X1)*(XP-X2) + (Y2-Y1)*(YP-Y2) + (Z2-Z1)*(ZP-Z2)

    CORE = 0.01 #to avoid singularities at vortex core
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
    #print(U,V,W)
    return [U,V,W]
#%%     CONTROL POINT DISTRIBUTION


b = 70     # wing span
n =  34  # Total number of points in the distribution
ncp = n-1   # control points - bound vortex

theta = np.linspace(0, np.pi, n)
cosine_values = np.cos(theta)

y_values = -b/2 * cosine_values + b/2    #uncheck for cosine distribution
distribution = "Cosine"

# y_values = np.linspace(0,b,n)           #uncheck for uniform distribution
# distribution = "Uniform"

# Define x and z coordinates (constant for wing)
x_values = np.zeros(ncp)
z_values = np.zeros(ncp)


# initialize the 3D coordinates for the bound vortex points.
coordinates_wing = np.zeros([ncp,3])

for i in range(ncp):
    coordinates_wing[i,1] = 0.5* (y_values[i]+y_values[i+1])  #bound vortex control points at the centre of the vortex filaments


#%%     OPERATING AND GEOMETRIC PARAMETERS

U0 = 1          #inflow vel
c = 1           #chord
rho = 1.225     #density


alpha = 5       #geometric aoa
alpha_rad = np.radians(alpha)

wake_steps = 20 #no of steps in wake

wake_distance = 20*b                #total wake distance to compute
wake_ds = wake_distance/wake_steps  #distance between steps, length of trailing vortex filaments

t = wake_ds/U0                      #timestep


#%%     INITIALIZE ARRAYS

# induced velocity matrix for unit trailing vortex circulation

# element (i,j) gives the induced velocity at control point i due to the trailing vortex j
u_matrix = np.zeros([ncp,n])
v_matrix = np.zeros([ncp,n])
w_matrix = np.zeros([ncp,n])

# Trailing vortices
gamma = (np.ones(n)).transpose() 
GammaNew_trail = gamma.copy()

# Bound vortices
alpha_eff_matrix = np.zeros(ncp)
CL_matrix = np.zeros(ncp)
Gamma_bound = np.zeros(ncp)
GammaNew_bound = Gamma_bound.copy()



#%%     TRAILING VORTEX POINTS

tvortex_pts = np.zeros([wake_steps, n, 3]) #trailing edge vortex points

# first set of trailing vortices are in the bound vortex line, ie quarter chord
tvortex_pts[0,:,0] = 0
tvortex_pts[:,:,1] = y_values
tvortex_pts[:,:,2] = 0

#For every trailing vortex,n, find the location of every wake_steps
for i in range(n): 
    for k in range(1,wake_steps):

        tvortex_pts[k,i,0] =  tvortex_pts[k-1,i,0] + U0*t
        
#%%     INDUCED VELOCITY MATRIX GENERATION
for i in range(ncp):
    for j in range(n):
        for k in range(wake_steps-1):
            
            # sign reversal of gamma is incorporated by interchanging the order of wake points. Right hand thumb rule.
            if j < n/2:
            
                # Compute induced velocity by each trailing filament
                Velocity = Biot_Savart(coordinates_wing[i], tvortex_pts[k+1,j,:],tvortex_pts[k,j,:]) 
                u = Velocity[0]
                v = Velocity[1]
                w = Velocity[2]
                
                #for a single trailing vortex and control point pair, all trailing filament contributions are added to a single element
                u_matrix[i][j] += u
                v_matrix[i][j] += v
                w_matrix[i][j] += w
            else:
                # Compute induced velocity by each trailing filament
                Velocity = Biot_Savart(coordinates_wing[i], tvortex_pts[k,j,:],tvortex_pts[k+1,j,:]) 
                u = Velocity[0]
                v = Velocity[1]
                w = Velocity[2]
                
                #for a single trailing vortex and control point pair, all trailing filament contributions are added to a single element
                u_matrix[i][j] += u
                v_matrix[i][j] += v
                w_matrix[i][j] += w
                
            

#%%     LIFTING LINE SOLVER

#max_itr = 1000
tol = 1e-12
itr = 0
error = 1

while np.all(error > tol):
    #print(itr)
    
    """
    u_matrix @ gamma gives the induced u velocity for actual gamma distribution
    
    """
    # Effective velocity at each control point for the actual trailing gamma distribution
    U_matrix = U0 + (u_matrix @ gamma)  
    V_matrix = v_matrix @ gamma
    W_matrix = w_matrix @ gamma #downwash
    
    V_eff_matrix = np.sqrt(U_matrix**2 + V_matrix**2 + W_matrix**2)
    
    
    alpha_i_matrix = np.arctan(W_matrix/ U_matrix) # induced AoA due to downwash, negative sign taken care in velocity direction
    alpha_eff_matrix = alpha_rad + alpha_i_matrix  # alpha_i_matrix sign depends on velocity direction, so add these two angles
    
    for j in range(ncp):
        
        # if itr == 0:
        #     dspan = tvortex_pts[0,j+1,1] - tvortex_pts[0,j,1]       # length of bound vortex filament
        dspan =1   
        alpha_eff = alpha_eff_matrix[j]                         # effective angle of attack
        V_eff = V_eff_matrix[j]                                 # effective velocity

        CL = compute_CL(alpha_eff)                              # CL at the bound filament

        CL_matrix[j] = CL                                       # CL distribution matrix
        
        
        Lift = 0.5*rho*(V_eff**2)*CL*c*(dspan)*np.cos(0.5*(alpha_i_matrix[j]))  #Lift - perpendicular to freestream
        
        gammaNew = Lift/(rho*V_eff)                             # Circulation from Lift values

        Gamma_bound[j] = gammaNew                               # new bound vortex from lift distribution
    
    
    for j in range(int(n/2)):                                   # compute trailing gamma from bound gamma
        if j == 0 :
            GammaNew_trail[j] = Gamma_bound[j]
            GammaNew_trail[n-1] = GammaNew_trail[j]

        else:
            GammaNew_trail[j] = Gamma_bound[j] - Gamma_bound[j-1]  
            GammaNew_trail[n-j-1] = GammaNew_trail[j]
        
    gamma = 0.75*gamma + 0.25*GammaNew_trail                    # trailing gamma estimate for next iteration
    
    
    
    print("iterations", itr)
    
    error = np.abs(gamma - GammaNew_trail)
    print("Max error:", np.max(error))
    itr += 1

    # if np.all(error < tol):
    #     print("solution converged")
    #     break
    
#%%     POST PROCESSING

fig = plt.figure(figsize=(12, 5))


ax = fig.add_subplot(121, projection='3d')


# Control points on the bound vortex
ax.scatter(coordinates_wing[:,0], coordinates_wing[:,1], coordinates_wing[:,2], color='r', alpha=1)

# Mesh for trailing vortex 
tvortex_x = tvortex_pts[:,:,0]
tvortex_y = tvortex_pts[:,:,1]
tvortex_z = tvortex_pts[:,:,2]
ax.plot_wireframe(tvortex_x, tvortex_y, tvortex_z)
ax.scatter(tvortex_pts[0,:,0],tvortex_pts[0,:,1],tvortex_pts[0,:,2])
ax.set_box_aspect([10, 3, 2]) # adjust the ratio of 3d plot
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.view_init(elev=30, azim=45)  # adjust  default viewing angle

# Subplot for CL distribution

plt.subplot(122)
plt.plot(coordinates_wing[:,1]/b, CL_matrix)

plt.title(f"CL distribution for {distribution} distribution with {int(n/2)} HS vortices")
plt.xlabel('Normalized Spanwise coordinate')
plt.ylabel('CL')
plt.ylim(0,np.max(CL_matrix)+0.1)
plt.grid()


plt.tight_layout()


plt.show()
