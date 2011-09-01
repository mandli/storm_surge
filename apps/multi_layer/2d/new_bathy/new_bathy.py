#!/usr/bin/env python
# encoding: utf-8

import sys

import numpy as np
import matplotlib.pyplot as plt

N = 10000
x = np.linspace(-200e3,500e3,N)
z = np.zeros(N)

x0 = 350e3
x1 = 450e3
x2 = 480e3
basin_depth = -3000
shelf_depth = -100
beach_slope = 0.05
eta = [0.0,-300]
h = 10.0

A = basin_depth - eta[1] + 0.5*h
B = shelf_depth - eta[1] - 0.5*h
eta_int = (A*x1 - B*x0) / (A-B)
if not (x0 < eta_int < x1):
    print "Intersection of internal interface does not hit continental slope."
    sys.exit(1)
    
shelf_slope =  A / (x0 - eta_int)

z += (x < x0) * basin_depth
z += (x0 < x) * (x < eta_int) * (shelf_slope * (x - x0) + basin_depth)
z += (eta_int < x) * (x < x1) * (shelf_slope * (x - x1) + shelf_depth)
z += (x1 < x) * (x < x2) * shelf_depth
z += (x2 < x) * ( beach_slope * (x - x2) + shelf_depth)

plt.hold(True)
layer = np.empty((N,2))
layer[:,0] = eta[0] * (eta[0] > z) + z * (eta[0] <= z)
layer[:,1] = eta[1] * (eta[1] > z) + z * (eta[1] <= z)
# plt.fill_between(x,z,eta[:,0],(eta[:,0] > z) * (eta[:,1] > z),color=(0.2,0.8,1.0))
plt.fill_between(x,z,layer[:,0],color=(0.2,0.8,1.0))
plt.fill_between(x,z,layer[:,1],color='b')
plt.plot(x,z,'k')
plt.plot(x,layer[:,0],'k')
plt.plot(x,layer[:,1],'k')
plt.axis([-200e3,500e3,basin_depth - 50,50])
plt.hold(False)
plt.show()