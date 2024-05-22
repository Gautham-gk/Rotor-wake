
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
    #print(U,V,W)
    return [U,V,W]


#%%   CONTROL POINT DISTRIBUTION


b = 10  # wing span
n = 10 # Total number of points in the distribution, no. of points that we are considering in each cross section


# Generate the cosine distribution for z values from 0 to pi
theta = np.linspace(0, np.pi, n)
cosine_values = np.cos(theta)


# Scale the cosine values to the wing span [-b/2, b/2]
#y_values = b/2 * cosine_values +b/2
y_values = np.linspace(0,b,n)

# Define x and z coordinates (constant for wing)
x_values = np.zeros(n)
z_values = np.zeros(n)


# Create the 3D coordinates. This is basically the coordinates of the quarter chord.
coordinates_wing = np.zeros([n,3])

coordinates_wing[:,1] = y_values


#%%  WAKE VORTICES DISTRIBUTION



U0 = 1 #inflow vel
c = 1 #chord

rho = 1.225


alpha = 5 #geometric aoa
alpha_rad = np.radians(alpha)

wake_steps = 20 #no of steps in wake, Number of y-z planes that we are considering

wake_distance = 20*b #total wake distance to compute
wake_ds = wake_distance/wake_steps #distance between steps

t = wake_ds/U0 #time taken between steps
# t is time for the wind particle to go through each wake




#initialize matrices
tvortex_pts = np.zeros([wake_steps, n, 3]) #trailing edge vortex points

# row 'n' represents the control points, and the column 'pair_num' represents the influence of trailing vortice on that point
u_matrix = np.zeros([n,n])
v_matrix = np.zeros([n,n])
w_matrix = np.zeros([n,n])


#The row sum would represent the total induced velocity at a point
u_matrix_sum = np.zeros(n)
v_matrix_sum = np.zeros(n)
w_matrix_sum = np.zeros(n)
V_eff_matrix_sum = np.zeros(n)

alpha_eff_matrix = np.zeros(n)
CL_matrix = np.zeros(n).transpose()

gamma = (np.ones(n)).transpose() #initialize gamma, we are assuming that gamma is 1 initially
#this gamma is at the control point, which is clockwise, so let it be, but change the w_matrix 
#remember two vortices of opposite signs.

# for i in range(n):
#     if i <n/2:
#         gamma[i] = -1
#     else:
#         gamma[i] = 1

GammaNew_bound = gamma.copy()
Gamma_bound = gamma.copy()
GammaNew_trail = gamma.copy()
Lift_list = np.zeros(n)
print("gamma initial:", gamma)

z = 0

for i in range(n):
    x = c #one chord length after the vortex filament #AKA- quarter chord after the TE
    for k in range(wake_steps):
        
        tvortex_pts[k,i,0] = -x
        
        x -= -U0*t
        
        tvortex_pts[k,i,1] = coordinates_wing[i,1]

for i in range(n):
    for j in range(n):
        for k in range(wake_steps-1):


            Velocity = Biot_Savart(coordinates_wing[i], tvortex_pts[k,j,:],tvortex_pts[k+1,j,:])
            u = Velocity[0]
            v = Velocity[1]
            w = Velocity[2]
            
            u_matrix[i][j] += u
            v_matrix[i][j] += v
            
            if j <n/2:
                if j < i:
                    w_matrix[i][j] +=  -np.abs(w)
                else:
                    w_matrix[i][j] +=  +np.abs(w)
            else:
                if j < i:
                    w_matrix[i][j] += +np.abs(w)
                else:
                    w_matrix[i][j] += -np.abs(w)

            

#%%       TRY TO PUT EVERYTHING UNDER ITERATING LOOP
#u_matrix has induced vel for unit circulation, hence multiplied with actual gamma
W_New =  w_matrix @ gamma
dspan_mat = []

max_itr = 100
tol = 1e-6
itr = 0
while itr < max_itr:
    print(itr)
    U_matrix = U0 + (u_matrix @ gamma)  #V_effective 
    V_matrix = v_matrix @ gamma
    W_matrix = w_matrix @ gamma
    V_eff_matrix = np.sqrt(U_matrix**2 + V_matrix**2 + W_matrix**2)
    alpha_i_matrix = np.arctan(W_matrix/ U_matrix)
    alpha_eff_matrix = alpha_rad + alpha_i_matrix
    
    for j in range(n):
        #dspan = coordinates_wing[j,1] - coordinates_wing[j-1,1]
        #dspan_mat.append(dspan)
        dspan = float(b)/float(n-1)
        #dspan = 1
        alpha_eff = 0.5*(alpha_eff_matrix[j] + alpha_eff_matrix[j-1])
        V_eff = 0.5*(V_eff_matrix[j] + V_eff_matrix[j-1])

        CL = compute_CL(alpha_eff)

        #print(CL)

        CL_matrix[j] = CL
        
        
        Lift = 0.5*rho*(V_eff**2)*CL*c*(dspan)*np.cos(0.5*(alpha_i_matrix[j]+ alpha_i_matrix[j-1]))
        gammaNew = Lift/(rho*V_eff)
        
        Lift_list[j] = Lift

        Gamma_bound[j] = gammaNew #new bound vortex from lift distribution
    #print("NEW GAMMA:",GammaNew_dist)
    
    for j in range(int(n/2)):
        if j == 0 :
            GammaNew_trail[j] = Gamma_bound[j]
            GammaNew_trail[n-1] = GammaNew_trail[j]

        else:
            GammaNew_trail[j] = Gamma_bound[j] - Gamma_bound[j-1]  
            GammaNew_trail[n-j-1] = GammaNew_trail[j]
        
    gamma = 0.75*gamma + 0.25*GammaNew_trail #this gamma is trailing vortices strength
    
    
    itr += 1
    #print("iterations", itr)
    error = np.abs(gamma - Gamma_bound) #this error value has to be changed
    print("error", error)
    
    print("effective alpha", alpha_eff_matrix)
    # print("gammabound", GammaNew_bound)
    if np.all(error < tol):
        print("solution converged")
        break



dpi_setting = 300
fig = plt.figure(dpi=dpi_setting)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_values, y_values, z_values, c=y_values, cmap='viridis')
ax.scatter(tvortex_pts[1:,:,0], tvortex_pts[1:,:,1], tvortex_pts[1:,:,2])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

plt.plot(y_values/b,CL_matrix)

