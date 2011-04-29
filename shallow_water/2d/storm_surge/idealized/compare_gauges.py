#!/usr/bin/env python

import os

import numpy as np
import matplotlib.pyplot as plt

# fig = plt.figure()
# ax = fig.add_subplot(111)
# 
# t = np.arange(0.0, 5.0, 0.01)
# s = np.cos(2*np.pi*t)
# line, = ax.plot(t, s, lw=2)
# 
# ax.annotate('local max', xy=(2, 1), xytext=(3, 1.5),
#             arrowprops=dict(facecolor='black', shrink=0.05),
#             )
# 
# ax.set_ylim(-2,2)
# plt.show()

from pyclaw.plotters.data import ClawPlotData

plotdata = ClawPlotData()

plt.figure(1)

plotdata.outdir = '_output'

plt.hold(True)
plt.grid(True)
gs = plotdata.getgauge(9)
plt.plot(gs.t,gs.q[:,3],'b',linewidth=2)
plt.annotate('Gauge 1',(100e3,2.0),(120e3,2.25),arrowprops={'width':1,'color':'k'})
gs = plotdata.getgauge(11)
plt.plot(gs.t,gs.q[:,3],'b--',linewidth=2)
plt.annotate('Gauge 2',(100e3,-2.0),(120e3,-2.25),arrowprops={'width':1,'color':'k'})
# gs = plotdata.getgauge(18)
# plt.plot(gs.t,gs.q[:,3],'b-.',linewidth=2)
# plt.annotate('Gauge 3',(50e3,-0.15),(15e3,-1.5),arrowprops={'width':1,'color':'k','shrink':0.1})    

plt.plot([0,40*3600],[0,0],'k')

    
plt.title("Gauge Data")
plt.xlabel('t (hours)')
plt.ylabel('m')

locs,labels = plt.xticks()
# import pdb; pdb.set_trace()
labels = np.trunc(locs/3600.0)
# locs = np.linspace(-12.0,40,52)
# labels = range(-12,41)
plt.xticks(locs,labels)
plt.axis([0,40*3600,-3,3])

plt.hold(False)

# plt.savefig('/Users/mandli/Documents/research/Papers/awr10/figures/ss_gauge.png')
plt.savefig('/Users/mandli/Desktop/ss_gauge.png')
plt.show()
    
