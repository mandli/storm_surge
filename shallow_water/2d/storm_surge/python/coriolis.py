#!/usr/bin/env python

import numpy as np
import matplotlib

params = {'axes.labelsize': 12,
          'text.fontsize': 12,
          'legend.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'text.usetex': True}
matplotlib.rcParams.update(params)

import matplotlib.pyplot as plt


# Parameters
N = 1000
R_earth = 6367500.0
omega = 7.2722052166430395e-05
theta_0 = 31.0
deg_range = 2.0
# theta_0 = 43.5
# deg_range = 43.0*2.0

y = np.linspace(-300e3,300e3,N)
theta = np.linspace(theta_0-deg_range/2.0,theta_0+deg_range/2.0,N)
theta_0 = theta_0

fig = plt.figure(1)
ax = fig.add_subplot(111)

exact = ax.plot(theta,2.0*omega*np.sin(theta * np.pi/180))
base = ax.plot(theta,np.ones(theta.shape)*2.0*omega*np.sin(theta_0 * np.pi/180))
y_theta = y/111e3 * np.pi/180 + theta_0 * np.pi/180
beta = ax.plot(theta,2.0*omega*(np.sin(theta_0 * np.pi/180) + (y_theta - theta_0 * np.pi/180) * np.cos(theta_0 * np.pi/180)))
# stommel = ax.plot(theta,1e-4 + 1e-11 * (y-0.0))

plt.legend( (exact,base,beta),  ('True f', 'Constant f', 'Beta-plane f'), loc=4)
# plt.legend( (exact,base,beta,stommel),  ('True f', 'Constant f', 'Beta-plane f','Stommel'), loc=4)


f = 2.0*omega*np.sin(theta_0 * np.pi/180)
L = np.linspace(50e3,500e3,100)
g = 9.81
r = 0.98

fig = plt.figure(2)
ax = fig.add_subplot(111)
rossby_plots = []
rossby_labels = ['$U = 2$',"$U = \sqrt{g 100}$","$U = \sqrt{g 200}$","$U = \sqrt{g 2000}$","$U = \sqrt{g 3000}$","$U = \sqrt{g 4000}$"]
# Could also use gravity wave speed = np.sqrt(g*4000.0)
for U in [2.0,np.sqrt(g*(1-r)*100.0),
          np.sqrt(g*(1-r)*200.0),
          np.sqrt(g*(1-r)*2000.0),
          np.sqrt(g*(1-r)*3000.0),
          np.sqrt(g*(1-r)*4000.0)]:
    Ro = U/(L*f)
    rossby_plots.append(ax.plot(L,Ro))
    print "Range(U = %s) = [%s,%s]" % (U,np.min(Ro),np.max(Ro))

ax.set_title("Rossby Number")
ax.set_xlabel('L (km)')
ax.set_ylabel('Ro')
locs,labels = plt.xticks()
labels = locs/1e3
plt.xticks(locs,labels)
ax.legend(rossby_plots,rossby_labels)
plt.savefig('rossby_numbers.pdf')

plt.show()