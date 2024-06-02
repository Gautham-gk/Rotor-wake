# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 01:20:14 2024

@author: gk220
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 10:08:46 2024

@author: Mathesh JK
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:00:34 2024

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

file = pd.read_excel(r"D:\TU Delft\Q3\Rotor Wake Aero\Assignment\polar DU95W180 (3).xlsx") #enter polar file location
# Extract columns for AoA, Cl, Cd, and Cm
AoA_values = file['Alfa'].tolist()
Cl_values = file['Cl'].tolist()
Cd_values = file['Cd'].tolist()
Cm_values = file['Cm'].tolist()

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

    CORE = 0.001 #to avoid singularities at vortex core
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
#%%     WAKE POINTS 

def wake_points_gen(a):
    
    u_rotor = U0*(1-a)
    u_far = U0*(1-2*a)    
    u_wake = np.linspace(u_rotor, u_far, wake_steps)
    time = np.linspace(0,wake_len/(0.5*(u_far+u_rotor)),wake_steps)
    
    #Modified to include the number of blades
    
    for num in range(N):
        for k in range(wake_steps):
            for i in range(n):
                nv[k,i,0,num] = +u_wake[k]*time[k]
                #tvortex[k,i,0] = 2*k
                
                nv[k,i,1,num] = bladey[i]*np.cos(omega*time[k] + num*b_offset)
                nv[k,i,2,num] = bladey[i]*np.sin(omega*time[k] + num*b_offset )
    return nv


#%%     INDUCED VELOCITY MATRIX
def induced_vel_mat(blade_coord,nv):
    for num in range(N):
        for i in range(ncp):
            for j in range(n):
                for k in range(wake_steps-1):
                    
                    # # Compute induced velocity by each trailing filament
                    # Velocity = Biot_Savart(blade_coordinates[i], tvortex_pts[k,j,:],tvortex_pts[k+1,j,:]) 
                    # u = Velocity[0]
                    # v = Velocity[1]
                    # w = Velocity[2]
                    
                    # sign reversal of gamma is incorporated by interchanging the order of wake points. Right hand thumb rule.
                    if j < n/2:
                    
                        # Compute induced velocity by each trailing filament
                        Velocity = Biot_Savart(blade_coord[i,:,num], nv[k+1,j,:,num],nv[k,j,:,num]) 
                        u = Velocity[0]
                        v = Velocity[1]
                        w = Velocity[2]
                        
                        #for a single trailing vortex and control point pair, all trailing filament contributions are added to a single element
                        u_matrix[i][j] += u
                        v_matrix[i][j] += v
                        w_matrix[i][j] += w
                    else:
                        # Compute induced velocity by each trailing filament
                        Velocity = Biot_Savart(blade_coord[i,:,num], nv[k,j,:,num],nv[k+1,j,:,num]) 
                        u = Velocity[0]
                        v = Velocity[1]
                        w = Velocity[2]
                        
                        #for a single trailing vortex and control point pair, all trailing filament contributions are added to a single element
                        u_matrix[i][j] += u
                        v_matrix[i][j] += v
                        w_matrix[i][j] += w
                    
    return u_matrix, v_matrix, w_matrix

#%%     DEFINE PARAMETERS

N=3 # No. of blades
b_offset= 2*np.pi/N
#b_offset is for finding x,y,z coordinates of different blades. FOr n=3, blades would be at 0,120 and 240 degrees.


#rotor parameters
root_mu = 0.2
tip_mu = 1
R = 50
c = 1
n = 16
ncp = n-1

#alpha = 5 #deg
#alpha_rad = np.deg2rad(alpha)


#mu_distribution = np.linspace(root_mu,tip_mu,n)
#This is for uniform

theta = np.linspace(np.pi, 0, n)
cosine_values = np.cos(theta)

mu_avg=0.5*(tip_mu+root_mu)
mu_distribution =0.5*(tip_mu-root_mu)*cosine_values + 0.5*(tip_mu+root_mu)





#operational parameters
U0 = 10
rho = 1.225
lam = 8

omega = lam*U0/R        #has to be calc from lambda

#wake parameters
wake_steps = 7000
wake_len = 5*R

#spanwise parameters

P = 2
Twist = [-14*(1-x) for x in mu_distribution]
B_distribution = [x + P for x in Twist]
chord_distribution = [3*(1-x)+1 for x in mu_distribution]




#%%     CONTROL POINT DISTRIBUTION


bladey = mu_distribution*R  #starting points of trailing vortices

blade_coord = np.zeros([ncp,3,N])
for num in range(N):
    for i in range(ncp):
         dist= 0.5*(bladey[i]+ bladey[i+1])
         blade_coord[i,1,num]=dist*np.cos(num*b_offset) #cos component is y coordinate
         blade_coord[i,2,num]=dist*np.sin(num*b_offset) #sin component is y coordinate



#%%     INITIALIZE ARRAYS

# induced velocity matrix for unit trailing vortex circulation

# element (i,j) gives the induced velocity at control point i due to the trailing vortex j
u_matrix = np.zeros([ncp,n])
v_matrix = np.zeros([ncp,n])
w_matrix = np.zeros([ncp,n])

a =0.25
    

# Trailing vortices
gamma = (np.ones(n)).transpose() 
GammaNew_trail = gamma.copy()

# Bound vortices
alpha_eff_matrix = np.zeros(ncp)
CL_matrix = np.zeros(ncp)
Gamma_bound = np.zeros(ncp)
GammaNew_bound = Gamma_bound.copy()


nv= np.zeros([wake_steps, n, 3,N])

#nv is basically t_vortex_pts, except that it a 4d matrix now, to account for the number of blades. 
# I changed the name of variable so that modifying the code would be easier.
#Similarly, I changed blade_coordinates to  blade_coord for same reason.



#new_tvortex = nv[wake_steps, n, 3,0] #trailing edge vortex points
#The above code started giving an error when modified, so I commented out.


nv = wake_points_gen(a)
u_matrix, v_matrix, w_matrix = induced_vel_mat(blade_coord, nv)

#%%     SOLVER

tol = 0.1
itr = 0
error = 1
a_prime = 0.01


itr = 0
max_itr = 200

while np.all(error > tol):
#while itr < max_itr:
    #print(itr)

    """
    u_matrix @ gamma gives the induced u velocity for actual gamma distribution
    
    """
    # Effective velocity at each control point for the actual trailing gamma distribution
    U_matrix = (u_matrix @ gamma)  
    V_matrix = v_matrix @ gamma
    W_matrix = w_matrix @ gamma   #downwash
    
    
    Um = -1* U_matrix
    #print(Um)
    Vaxial = U0 + Um
    Vazim=0
    # for num in range(N):
    #     nazim=np.cross([1,0,0],blade_coord[1,:,num])
        
        # Normalize the vector to get the unit normal vector
        # norm = np.linalg.norm(nazim)
        # unit_nazim = nazim / norm
        
        
    Vazim = W_matrix + omega*blade_coord[:,1,0] 
    # Dont know how to take into account the different blades here, but still getting good results.
    
    V_eff_matrix = np.sqrt(Vaxial**2 +  Vazim**2)
    
    
    phi_matrix = np.arctan(Vaxial/ Vazim) # induced AoA due to downwash, negative sign taken care in velocity direction
    alpha_eff_matrix = np.rad2deg( phi_matrix ) # alpha_i_matrix sign depends on velocity direction, so add these two angles
    
    for j in range(ncp):
        
        mu = blade_coord[j,1,0]/R
        A = np.pi*((R*mu_distribution[j+1])**2 - (R*mu_distribution[j])**2) 
        B = np.interp(mu, mu_distribution, B_distribution)
        c = np.interp(mu, mu_distribution, chord_distribution)
        
        #print("B", B)
        
        alpha_eff = alpha_eff_matrix[j]  + B                       # effective angle of attack
        alpha_eff_matrix[j] = alpha_eff
        phi = phi_matrix[j]
        V_eff = V_eff_matrix[j]                                 # effective velocity

        Cl = np.interp(alpha_eff, AoA_values, Cl_values)        # Cl at the bound filament
        Cd = np.interp(alpha_eff, AoA_values, Cd_values)

        CL_matrix[j] = Cl                                      # CL distribution matrix
        
        
        Lift = 0.5*rho*(V_eff**2)*Cl*c                          #Lift 
        Drag = 0.5*rho*(V_eff**2)*Cd*c                          # Drag
        
        Faxial = Lift * np.cos(phi) + Drag*np.sin(phi)
        Fazim = Lift * np.sin(phi) - Drag*np.cos(phi)
        
        Thrust = Faxial*(R*(mu_distribution[j+1]-mu_distribution[j]))
        Torque = Fazim*(R*(mu_distribution[j+1]-mu_distribution[j]))
        
        CT = (Thrust*N)/(0.5*rho*(U0**2)*A)
        
        
        
        gammaNew = Lift/(rho*V_eff)                             # Circulation from Lift values

        Gamma_bound[j] = gammaNew                               # new bound vortex from lift distribution
    
    
    for j in range(int(n/2)):                                   # compute trailing gamma from bound gamma
        if j == 0 :
            GammaNew_trail[j] = Gamma_bound[j]
            GammaNew_trail[n-1] = GammaNew_trail[j]

        else:
            GammaNew_trail[j] = Gamma_bound[j] - Gamma_bound[j-1]  
            GammaNew_trail[n-j-1] = GammaNew_trail[j]
    
    m = 0.6
    gamma = m*gamma + (1-m)*GammaNew_trail                    # trailing gamma estimate for next iteration
    #itr += 1
    
    
    print("iterations", itr)
    
    error = np.abs(gamma - GammaNew_trail)
    print("Max error:", np.max(error))
    print("min error",np.min(error))
    itr += 1

    # if np.all(error < tol):
    #     print("solution converged")
    #     break

#%%     POST PROCESSING
fig = plt.figure(figsize=(12, 5))


ax = fig.add_subplot(121, projection='3d')


# I am plotting only the vortex from one blade
# Control points on the bound vortex
ax.scatter(blade_coord[:,0,0], blade_coord[:,0,1], blade_coord[:,0,2], color='r', alpha=1)

# Mesh for trailing vortex 


tvortex_x = nv[:,:,0,num]
tvortex_y = nv[:,:,1,num]
tvortex_z = nv[:,:,2,num]
ax.scatter(nv[0,:,0,num],nv[0,:,1,num],nv[0,:,2,num])
ax.plot_wireframe(tvortex_x, tvortex_y, tvortex_z)

#ax.set_box_aspect([10, 3, 2]) # adjust the ratio of 3d plot
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.view_init(elev=30, azim=45)  # adjust  default viewing angle

# Subplot for CL distribution

plt.subplot(122)
plt.plot(blade_coord[:,1,0]/R, np.rad2deg(phi_matrix),marker='o',label='phi')
plt.plot(blade_coord[:,1,0]/R, alpha_eff_matrix ,marker='o',label='AOA')

plt.title(f"CL distribution for constant distribution with {int(n/2)} HS vortices")
plt.xlabel('Normalized Spanwise coordinate')
plt.ylabel('CL')
#plt.ylim(0,np.max(CL_matrix)+0.1)
plt.grid()

plt.legend()
plt.tight_layout()


plt.show()