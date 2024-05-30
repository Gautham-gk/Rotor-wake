# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:42:10 2024

@author: gk220
"""

import numpy as np
import math
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d

U0 =10
R=50 #Rotor radius
a_w=0.2 #Assumed
omega= 0.2 #assumed
U_w=U0*(1-a_w)
U_far = U0*(1-2*a_w)


wake_steps = 20 #no of steps in wake, Number of y-z planes that we are considering
b = 10  # wing span
n = 10
#n is number of control points

wake_distance = 20*b #total wake distance to compute
wake_ds = wake_distance/wake_steps #distance between steps
t_incr = wake_ds/U0 #time taken between steps

# t is time for the wind particle to go through each wake

U_wake = np.linspace(U_w,U_far,wake_steps)

Radius_wake = np.linspace(R,np.sqrt(U_w/U_far*R**2),wake_steps)
#Wtf is this?



root_mu=0.2
tip_mu=1

y_values = np.linspace(root_mu,tip_mu,n)
coordinates_wing = np.zeros([n,3])

coordinates_wing[:,1] = y_values
# Convert the NumPy array to a list
coordinates_wing = coordinates_wing.tolist()
print(coordinates_wing)

all_coordinates=[]
all_coordinates.append(coordinates_wing)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

M=1000

# Plot initial coordinates with same color for points with same j value
colors = plt.cm.viridis(range(n))
for j, coords in enumerate(coordinates_wing):
    x, y, z = coords
    ax.scatter(x, y, z, color=colors[j], marker='o')
t=0
for i in range(M):
    add_coord = []
    for j in range(n):
        #print("i and j", i, j)
        r=y_values[j]
        #print("i j and y_values",i,j)
        x_last = all_coordinates[i][j][0]
        y_last = all_coordinates[i][j][1]
        z_last = all_coordinates[i][j][2]
        #print("x_last,y_last,z_last", x_last, y_last, z_last)
        x_new = x_last+ (U_w*t)
        y_new =  r * math.sin(omega * t)
        z_new =  r * math.cos(omega * t)
        
        new_point = [x_new, y_new, z_new]  # Creating the new point as a list
        add_coord.append(new_point)
        
        # Define color based on j value
        color = plt.cm.viridis(j / n)
        
        # Plot new point
        ax.plot([x_last, x_new], [y_last, y_new], [z_last, z_new], color=color, marker='o', linestyle='-')
    
    t+=t_incr
    all_coordinates.append(add_coord)
    #print(all_coordinates)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
