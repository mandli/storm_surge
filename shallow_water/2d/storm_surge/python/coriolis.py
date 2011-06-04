#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 1000
R_earth = 6367500.0
omega = 7.2722052166430395e-05
theta_0 = 30.0
deg_range = 4.0

y = np.linspace(-300e3,300e3,N)
theta = np.linspace(theta_0-deg_range/2.0,theta_0+deg_range/2.0,N) * np.pi / 180
theta_0 = theta_0 * np.pi/180

fig = plt.figure(1)
ax = fig.add_subplot(111)

exact = ax.plot(theta,2.0*omega*np.sin(theta))
base = ax.plot(theta,np.ones(theta.shape)*2.0*omega*np.sin(theta_0))
y = y/111e3 * np.pi/180 + theta_0
beta = ax.plot(theta,2.0*omega*(np.sin(theta_0) + (y - theta_0) * np.cos(theta_0)))

plt.legend( (exact, base, beta),  ('True f', 'Constant f', 'Beta-plane f'), loc=4)
plt.show()