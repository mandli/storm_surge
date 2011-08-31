#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Hurricane parameters
A = 23.0           # Hurricane model fit parameter
B = 1.5     
Pn = 1005.0        # Nominal atmospheric pressure     
Pc = 950.0         # Pressure in the eye of the hurricane    
rho_air = 1.15     # Density of air

N = 1000
ramp_up_time = 60*60*12.0
r = np.linspace(0,50e3,N)
time = np.linspace(-ramp_up_time,60*60*2.0,9)
p = np.empty((N,len(time)))
for (i,t) in enumerate(time):
    p[:,i] = Pc + (Pn - Pc) * np.exp(-A/(r/1e3)**B)
    if t < 0.0:
        p[:,i] = Pn + (p[:,i]-Pn) * np.exp(-(t/(ramp_up_time*0.45))**2)

#
for (i,t) in enumerate(time):
    plt.subplot(3,3,i+1)
    plt.plot(r,p[:,i])    
    plt.axis([0,50e3,950,1010])
    plt.title('t=%s' % t)
plt.show()