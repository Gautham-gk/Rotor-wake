# -*- coding: utf-8 -*-
"""
Created on Thu May 30 16:47:18 2024

@author: Mathesh JK
"""

import numpy as np
import matplotlib.pyplot as plt

u = 1
omega = 1
a = 0.1
R = 50
r = 10

u_rotor = u*(1-a)
u_far = u*(1-2*a)

wake_len = 2*R
wake_steps = 5000
u_wake = np.linspace(u_rotor, u_far, wake_steps)
time = np.linspace(0,wake_len/(0.5*(u_far+u_rotor)),wake_steps)
n = 6



bladey = np.linspace(r,R,n)

tvortex = np.zeros([wake_steps,n,3])

for k in range(wake_steps):
    t = 0
    for i in range(n):
        #tvortex[k,i,0] = -u_wake[k]*time[k]
        tvortex[k,i,0] = 2*k
        
        tvortex[k,i,1] = bladey[i]*np.cos(omega*time[k])
        tvortex[k,i,2] = bladey[i]*np.sin(omega*time[k])
        
        

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(tvortex[:,:,0], tvortex[:,:,1], tvortex[:,:,2])
ax.set_box_aspect([10, 3, 2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')