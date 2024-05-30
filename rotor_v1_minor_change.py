# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:00:34 2024

@author: Mathesh JK
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#%%     AIRFOIL DATA
file = pd.read_excel(r"C:\Users\sunny\Documents\TU Delft AWE Course Stuffs\Q3 Courses and Books\AE4135 Rotor-Wake Aerodynamics\polar DU95W180 (3).xlsx") #enter polar file location
#C:\Users\sunny\Documents\TU Delft AWE Course Stuffs\Q3 Courses and Books\AE4135 Rotor-Wake Aerodynamics\BEM_turbine\polar_DU95W180.xlsx"
# Extract columns for AoA, Cl, Cd, and Cm
AoA_values = file['Alfa'].tolist()
Cl_values = file['Cl'].tolist()
Cd_values = file['Cd'].tolist()
Cm_values = file['Cm'].tolist()




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
#%%     WAKE POINTS 

def wake_points_gen(a):
    
    u_rotor = U0*(1-a)
    u_far = U0*(1-2*a)    
    u_wake = np.linspace(u_rotor, u_far, wake_steps)
    time = np.linspace(0,wake_len/(0.5*(u_far+u_rotor)),wake_steps)
    
    for k in range(wake_steps):
        for i in range(n):
            tvortex_pts[k,i,0] = -u_wake[k]*time[k]
            #tvortex[k,i,0] = 2*k
            
            tvortex_pts[k,i,1] = bladey[i]*np.cos(omega*time[k])
            tvortex_pts[k,i,2] = bladey[i]*np.sin(omega*time[k])
    return tvortex_pts

#%%     INDUCED VELOCITY MATRIX
def induced_vel_mat(blade_coordinates,tvortex_pts):
    for i in range(ncp):
        for j in range(n):
            for k in range(wake_steps-1):
                
                # Compute induced velocity by each trailing filament
                Velocity = Biot_Savart(blade_coordinates[i], tvortex_pts[k,j,:],tvortex_pts[k+1,j,:]) 
                u = Velocity[0]
                v = Velocity[1]
                w = Velocity[2]
                
                #for a single trailing vortex and control point pair, all trailing filament contributions are added to a single element
                u_matrix[i][j] += u
                v_matrix[i][j] += v
                
                if j <n/2:  # trailing vortices that have CW circulation
                    if j ==i:
                        
                        w_matrix[i][j] +=  -np.abs(w)
                    elif i <j:
                        w_matrix[i][j] +=  +np.abs(w)
                    else:
                        w_matrix[i][j] +=  -np.abs(w)
                        
                else:       # trailing vortices that have CCW circulation
                    if i == j:
                        w_matrix[i][j] +=  +np.abs(w)
                    elif i <j:
                        w_matrix[i][j] +=  -np.abs(w)
                    else:
                        w_matrix[i][j] +=  +np.abs(w)
    return u_matrix, v_matrix, w_matrix
#%%     DEFINE PARAMETERS

#rotor parameters
root_mu = 0.2
tip_mu = 1
R = 50
c = 1
n = 16
ncp = n-1
alpha = 5 #deg
alpha_rad = np.deg2rad(alpha)

#operational parameters
U0 = 10
rho = 1.225
lam=6 

omega = lam*U0/R      #has to be calc from lambda
print(omega)
#wake parameters
wake_steps = 5000
wake_len = 2*R




#%%     CONTROL POINT DISTRIBUTION
mu_distribution = np.linspace(root_mu,tip_mu,n)


bladey = mu_distribution*R  #starting points of trailing vortices

blade_coordinates = np.zeros([ncp,3])

for i in range(ncp):
    blade_coordinates[i,1] = 0.5*(bladey[i]+ bladey[i+1])


#spanwise parameters
P = 2
Twist = [-14*(1 - x) for x in mu_distribution] #degrees, blade twist
B_distribution = [x + P for x in Twist] #degrees, total pitch
chord_distribution = [3*(1-x) + 1 for x in mu_distribution] #m, chord distribution


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

tvortex_pts = np.zeros([wake_steps, n, 3]) #trailing edge vortex points

#%%     SOLVER

tol = 1e-3
itr = 0
error = 1
a = 0.2
a_prime = 0.2

while np.all(error > tol):
    #print(itr)
    
    tvortex_pts = wake_points_gen(a)
    u_matrix, v_matrix, w_matrix = induced_vel_mat(blade_coordinates, tvortex_pts)
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

#%% V azim 
# V_azim = [Uinf + u v w].n_azim +omega*r. n_azim is the cross product of omega and cooordinate_wing. So, it is i X j= k, taking unite normal vector, k_cap=n_azim
# therefore, V_azim = omega*r + w

Faxim_val=[]
Fazim_val=[]
CT_val=[]
Alfa_val=[]
phi_val=[]
CP_val=[]

for i in range(ncp):
    
    mu_mean= blade_coordinates[i]/R
    B = np.interp(mu_mean, mu_distribution, B_distribution) #deg, mean beta
    c = np.interp(mu_mean, mu_distribution, chord_distribution) #m, mean chord
    
    
    Ur= omega*blade_coordinates[i] + W_matrix[i]
    Ua= U_matrix[i]
    
    
    W = np.sqrt(Ur**2 + Ua**2) #m/s, inflow speed
    phi = np.arctan(Ur/Ua) #rad, inflow angle
    phi_deg =np.degrees(phi) #deg, inflow angle
    phi_val.append(phi_deg)
    Alfa = phi_deg + B #deg, angle of attack
    Alfa_val.append(Alfa)
    
            
    #extract Cl and Cd from the excel file data
    Cl = np.interp(Alfa, AoA_values, Cl_values)
    Cd = np.interp(Alfa, AoA_values, Cd_values)
            
    Lift = 0.5*rho*(W**2)*c*Cl #N, 3D lift force by the blade element
    Drag = 0.5*rho*(W**2)*c*Cd #N, 3D drag force by the blade element
            
    Faxial = (Lift*np.cos(phi) + Drag*np.sin(phi)) #N, axial force
    Fazim = (Lift*np.sin(phi) - Drag*np.cos(phi))   #N, azimuthal force
    Faxim_val.append(Faxial)
    Fazim_val.append(Fazim)
    
    if i < ncp - 1:
        Thrust = Faxial * (R * (mu_distribution[i+1] - mu_distribution[i]))
        Torque = Fazim * (R * (mu_distribution[i+1] - mu_distribution[i])) * R * mu_mean
        
        A = np.pi * ((bladey[i+1])**2 - (bladey[i])**2)
        CT = (Thrust) / (0.5 * rho * (U0**2) * A)  # Thrust coefficient from momentum theory
        CT_val.append(CT)
        
        Power = Torque * omega
        CP = Power / (0.5 * rho * (U0**3) * A)
        CP_val.append(CP)
        
    # Thrust = Faxial*(R*(mu_distribution[i+1]-mu_distribution[i]))
    # Torque = Fazim*(R*(mu_distribution[i+1]-mu_distribution[i]))*R*mu_mean
            
    # A = np.pi*((blade_coordinates[i+1])**2 - (blade_coordinates[i])**2)
    # CT = (Thrust)/(0.5*rho*(U0**2)*A) #Thrust coefficient from momentum theory
    # CT_val.append(CT)
    
    # Power = Torque*omega
    # CP = Power / (0.5*rho*(U0**3)*A)
    # CP_val.append(CP)
    


# In CT,CP formula,add N to account for number of blades

#%%     POST PROCESSING
fig = plt.figure(figsize=(12, 5))


ax = fig.add_subplot(121, projection='3d')


# Control points on the bound vortex
ax.scatter(blade_coordinates[:,0], blade_coordinates[:,1], blade_coordinates[:,2], color='r', alpha=1)

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
plt.plot(blade_coordinates[:,1]/R, CL_matrix)

plt.title(f"CL distribution for constant distribution with {int(n/2)} HS vortices")
plt.xlabel('Normalized Spanwise coordinate')
plt.ylabel('CL')
plt.ylim(0,np.max(CL_matrix)+0.1)
plt.grid()


plt.tight_layout()


plt.show()