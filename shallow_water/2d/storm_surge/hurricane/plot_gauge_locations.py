#!/usr/bin/env python
r"""Simple script that plots the gauge locations based on the data files
./amr2ez.data
./setgauges.data
"""

import numpy as np
import matplotlib.pyplot as plt

import pyclaw.data

run_data = pyclaw.data.Data('./amr2ez.data')

# Extract domain parameters
mx = run_data.mx
my = run_data.my
x = np.linspace(run_data.xlower,run_data.xupper,mx)
y = np.linspace(run_data.ylower,run_data.yupper,my)

# Read in gauge data
gauge_file = open('./setgauges.data','r')
# Read through header
for i in xrange(6):
    gauge_file.readline()
num_gauges = int(gauge_file.readline().split()[0])
gauges = np.empty((num_gauges,4))
for i in xrange(num_gauges):
    data = gauge_file.readline().split()
    gauges[i,:] = [data[1],data[2],data[3],data[4]]
    
fig = plt.figure(1)
ax = fig.add_subplot(111)
# Plot gauges
for i in xrange(num_gauges):
    ax.plot(gauges[i,0],gauges[i,1],'ko',markersize=5)
    ax.text(gauges[i,0]-35e3,gauges[i,1]-10e3,str(i).zfill(2),fontsize=15)
# Plot bathy locations
for bathy_ref in [350e3,450e3,480e3]:
    ax.plot(bathy_ref*np.ones(y.shape),y,'k--')
    
ax.set_title("Gauge Locations")
ax.set_xbound((run_data.xlower,run_data.xupper))
ax.set_ybound((run_data.ylower,run_data.yupper))
ax.set_xticks([-200e3,-100e3,0,100e3,200e3,300e3,400e3,500e3])
ax.set_yticks([-200e3,-100e3,0,100e3,200e3])
ax.set_xticklabels([-200,-100,0,100,200,300,400,500])
ax.set_yticklabels([-200,-100,0,100,200])
ax.set_xlabel('km')
ax.set_ylabel('km')
plt.grid(True)

plt.savefig("gauge_locations.pdf")
plt.show()