#!/usr/bin/env python
r"""Plots of new hurricane wind and pressure profiles"""

import sys

import numpy as np
import matplotlib.pyplot as plt

# Hurricane parameters
omega = 2.0 * np.pi / (60 **2 * 24) # Angular speed of Earth
rho_air = 1.15     # Density of air (kg/m^3)
V_max = 20.0       # Radius of maximum winds (m/s)
R_vmax = 7e3       # Radius of maximum winds (m)
p_n = 1005         # Nominal pressure (mb)
p_c = 950          # Central pressure (mb)
phi = 30.0         # Central latitude
f = 2.0 * omega * np.sin(phi * np.pi / 180.0)            # Coriolis parameter
# R_vmax = 4.6e1 * np.exp(-1.55e-2 * V_max + 1.69e-2 * phi * np.pi / 180.0) # Alternative derivation

# Fit parameters
B = V_max**2 * rho_air * np.exp(1.0) / (p_n - p_c) * 1.1e-2
B = 1.5
A = (R_vmax*1e-3)**(B)
# A = 23.0
print A,B

def wind_profile(r):
    x = r*1e-3
    return (np.sqrt( 1e1**2 * A * B * (p_n - p_c) * np.exp(-A/x**B) / (rho_air * x**B) 
        + x**2 * f**2 / 4.0) - x * f / 2.0)

def pressure_profile(r):
    x = r*1e-3
    return p_c + (p_n - p_c) * np.exp(-A/x**B)

N = 1000
r = np.linspace(0,500e3,N)
w = wind_profile(r)
p = pressure_profile(r)


fig = plt.figure(1)
ax = fig.add_subplot(211)
ax.plot(r,w)
ax.set_title('Wind Profile')
ax.set_xlim([0,50e3])
loc,label = plt.xticks()
label = loc/1.e3
plt.xticks(loc,label)
plt.xlabel('km')
plt.xticks(loc,label)
plt.xlabel('km')

ax = fig.add_subplot(212)
ax.plot(r,p)
ax.set_title('Pressure Profile')
ax.set_xlim([0,50e3])
loc,label = plt.xticks()
label = loc/1.e3
plt.xticks(loc,label)
plt.xlabel('km')

plt.show()